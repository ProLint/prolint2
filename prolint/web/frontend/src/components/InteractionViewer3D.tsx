/** 3D viewer for protein-ligand interactions using Molstar. */

import { useEffect, useRef, memo } from 'react';
import { Box, Typography, CircularProgress, Alert } from '@mui/material';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Color } from 'molstar/lib/mol-util/color';
import 'molstar/lib/mol-plugin-ui/skin/light.scss';
import { useMolstarPlugin, hexToMolColor } from '../hooks/useMolstarPlugin';
import { colors } from '../theme/visualizationTheme';

interface InteractionViewer3DProps {
  pdbData: string | null;
  queryResidue: number;
  databaseResidue?: number | null;
  loading?: boolean;
  error?: string | null;
}

function InteractionViewer3DComponent({
  pdbData,
  queryResidue,
  databaseResidue,
  loading,
  error,
}: InteractionViewer3DProps) {
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

  // Track previous data to avoid unnecessary reloads
  const prevDataRef = useRef<string | null>(null);

  // Create a key from the data that affects rendering - use content sample, not just length
  const dataKey = pdbData ? `${pdbData.slice(0, 200)}-${queryResidue}-${databaseResidue}` : null;

  // Load structure when data changes
  useEffect(() => {
    if (!isReady || !plugin || !pdbData || dataKey === prevDataRef.current) {
      return;
    }

    prevDataRef.current = dataKey;
    const currentLoadId = ++loadIdRef.current;

    const loadStructure = async () => {
      setIsLoading(true);

      try {
        const result = await loadPdbData(pdbData);
        if (!result || loadIdRef.current !== currentLoadId) return;

        const { structure } = result;

        // Add protein backbone as semi-transparent cartoon
        const polymer = await plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');
        if (polymer && loadIdRef.current === currentLoadId) {
          await plugin.builders.structure.representation.addRepresentation(polymer, {
            type: 'cartoon',
            color: 'uniform',
            colorParams: { value: Color(hexToMolColor(colors.neutral[400])) },
            typeParams: { alpha: 0.7 },
          });
        }

        if (loadIdRef.current !== currentLoadId) return;

        // Highlight query residue
        const queryComponent = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure,
          MS.struct.generator.atomGroups({
            'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), queryResidue]),
            'group-by': MS.struct.atomProperty.macromolecular.residueKey(),
          }),
          `query-${queryResidue}`
        );

        if (queryComponent && loadIdRef.current === currentLoadId) {
          // Cartoon + ball-and-stick for query residue
          await Promise.all([
            plugin.builders.structure.representation.addRepresentation(queryComponent, {
              type: 'cartoon',
              color: 'uniform',
              colorParams: { value: Color(hexToMolColor(colors.data.query)) },
            }),
            plugin.builders.structure.representation.addRepresentation(queryComponent, {
              type: 'ball-and-stick',
              color: 'uniform',
              colorParams: { value: Color(hexToMolColor(colors.data.query)) },
            }),
          ]);
        }

        if (loadIdRef.current !== currentLoadId) return;

        // Highlight database residue if provided
        if (databaseResidue) {
          const dbComponent = await plugin.builders.structure.tryCreateComponentFromExpression(
            structure,
            MS.struct.generator.atomGroups({
              'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), databaseResidue]),
              'group-by': MS.struct.atomProperty.macromolecular.residueKey(),
            }),
            `database-${databaseResidue}`
          );

          if (dbComponent && loadIdRef.current === currentLoadId) {
            await plugin.builders.structure.representation.addRepresentation(dbComponent, {
              type: 'ball-and-stick',
              color: 'uniform',
              colorParams: { value: Color(hexToMolColor(colors.data.database)) },
              typeParams: { sizeFactor: 0.3 },
            });
          }
        }

        if (loadIdRef.current !== currentLoadId) return;

        // Add ligands
        const ligand = await plugin.builders.structure.tryCreateComponentStatic(structure, 'ligand');
        if (ligand && loadIdRef.current === currentLoadId) {
          await plugin.builders.structure.representation.addRepresentation(ligand, {
            type: 'ball-and-stick',
            color: 'uniform',
            colorParams: { value: Color(hexToMolColor(colors.data.databaseLight)) },
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
  }, [dataKey, isReady, plugin, pdbData, queryResidue, databaseResidue, loadPdbData, setIsLoading, setViewerError, loadIdRef]);

  const showOverlay = loading || isLoading || error || viewerError || !pdbData;

  return (
    <Box>
      <Box sx={{ position: 'relative' }}>
        <div
          ref={containerRef}
          style={{
            width: '100%',
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
              borderRadius: '8px',
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
        <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mt: 1, mb: 0, textAlign: 'center' }}>
          Drag to rotate | Scroll to zoom | Right-drag to pan
        </Typography>
      )}
    </Box>
  );
}

export default memo(InteractionViewer3DComponent);
