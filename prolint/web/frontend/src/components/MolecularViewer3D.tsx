/** 3D molecular viewer with B-factor coloring using Molstar. */

import { useEffect, useRef, memo } from 'react';
import { Box, Typography, CircularProgress, Alert } from '@mui/material';
import 'molstar/lib/mol-plugin-ui/skin/light.scss';
import { useMolstarPlugin } from '../hooks/useMolstarPlugin';
import { colors } from '../theme/visualizationTheme';

interface MolecularViewer3DProps {
  pdbData: string | null;
  metric: string;
  loading?: boolean;
  error?: string | null;
}

function MolecularViewer3DComponent({ pdbData, metric: _metric, loading, error }: MolecularViewer3DProps) {
  const {
    containerRef,
    plugin,
    isReady,
    isLoading,
    error: viewerError,
    loadIdRef,
    setIsLoading,
    setError: setViewerError,
    loadPdbData,
  } = useMolstarPlugin();

  // Track previous pdbData to avoid unnecessary reloads
  const prevPdbDataRef = useRef<string | null>(null);

  // Load structure when pdbData changes
  useEffect(() => {
    // Skip if not ready, no data, or same data
    if (!isReady || !plugin || !pdbData || pdbData === prevPdbDataRef.current) {
      return;
    }

    prevPdbDataRef.current = pdbData;
    const currentLoadId = ++loadIdRef.current;

    const loadStructure = async () => {
      setIsLoading(true);

      try {
        const result = await loadPdbData(pdbData);
        if (!result || loadIdRef.current !== currentLoadId) return;

        const { structure } = result;

        // Fast cartoon representation with B-factor coloring
        const polymer = await plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');
        if (polymer && loadIdRef.current === currentLoadId) {
          await plugin.builders.structure.representation.addRepresentation(polymer, {
            type: 'cartoon',
            color: 'uncertainty',
            colorParams: {
              palette: {
                name: 'colorbrewer-rdylbu',
                params: { list: 'RdYlBu' },
              },
            },
          });
        }

        // Reset camera
        if (loadIdRef.current === currentLoadId) {
          plugin.managers.camera.reset();
          setViewerError(null);
        }
      } catch (err: unknown) {
        if (loadIdRef.current === currentLoadId) {
          const errorMessage = err instanceof Error ? err.message : String(err);
          setViewerError('Failed to load structure: ' + errorMessage);
        }
      } finally {
        if (loadIdRef.current === currentLoadId) {
          setIsLoading(false);
        }
      }
    };

    loadStructure();
  }, [pdbData, isReady, plugin, loadPdbData, setIsLoading, setViewerError, loadIdRef]);

  const showOverlay = loading || isLoading || error || viewerError || !pdbData;

  return (
    <Box>
      <Box sx={{ position: 'relative' }}>
        <div
          ref={containerRef}
          style={{
            width: '92%',
            aspectRatio: '1 / 1',
            minHeight: 350,
            border: `1px solid ${colors.border.default}`,
            borderRadius: '8px',
            overflow: 'hidden',
          }}
        />

        {showOverlay && (
          <Box
            sx={{
              position: 'absolute',
              top: 0,
              left: 0,
              right: 0,
              bottom: 0,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              backgroundColor: 'rgba(255, 255, 255, 0.9)',
              borderRadius: '4px',
              zIndex: 10,
            }}
          >
            {(loading || isLoading) && <CircularProgress />}
            {error && !loading && <Alert severity="error">{error}</Alert>}
            {viewerError && !loading && !error && <Alert severity="error">{viewerError}</Alert>}
            {!pdbData && !loading && !isLoading && !error && !viewerError && (
              <Alert severity="info">No structure data available</Alert>
            )}
          </Box>
        )}
      </Box>

      {pdbData && !error && !viewerError && (
        <Typography variant="caption" color="text.secondary" sx={{ mt: 1, mb: 0, display: 'block', textAlign: 'center' }}>
          Drag to rotate | Scroll to zoom | Right-drag to pan
        </Typography>
      )}
    </Box>
  );
}

export default memo(MolecularViewer3DComponent);
