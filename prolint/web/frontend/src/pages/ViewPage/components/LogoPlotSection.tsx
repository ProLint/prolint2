/** Logo plot visualization for residue selection. */

import { memo } from 'react';
import { Box, Typography, Fade } from '@mui/material';
import { colors, colorScales } from '../../../theme/visualizationTheme';
import ResidueBox from './ResidueBox';
import { LogoPlotSkeleton } from './Skeletons';
import type { LogoPlotData } from '../types';

interface LogoPlotSectionProps {
  logoPlotData: LogoPlotData | null;
  loading: boolean;
  selectedMetric: string;
  selectedResidues: number[];
  logoPlotNormalization: {
    minVal: number;
    maxVal: number;
    normalizeValue: (val: number) => number;
  };
  residuesByChain: Record<string, LogoPlotData['residues']>;
  isDraggingRef: React.MutableRefObject<boolean>;
  visualSelectionRef: React.MutableRefObject<Set<number>>;
  onMouseDown: (resid: number) => void;
  onMouseEnter: (resid: number) => void;
  onMouseUp: (resid: number) => void;
  onMouseLeave: () => void;
}

const LogoPlotSection = memo(function LogoPlotSection({
  logoPlotData,
  loading,
  selectedMetric,
  selectedResidues,
  logoPlotNormalization,
  residuesByChain,
  isDraggingRef,
  visualSelectionRef,
  onMouseDown,
  onMouseEnter,
  onMouseUp,
  onMouseLeave,
}: LogoPlotSectionProps) {
  if (loading) return <LogoPlotSkeleton />;

  if (!logoPlotData || logoPlotData.residues.length === 0) {
    return <Typography color="text.secondary">No data available</Typography>;
  }

  const { normalizeValue, minVal, maxVal } = logoPlotNormalization;
  const chainIds = Object.keys(residuesByChain).sort();
  const metricLabel = selectedMetric === 'occupancy' ? 'Occupancy' : selectedMetric === 'mean' ? 'Mean Duration' : 'Max Duration';
  const currentSelection = isDraggingRef.current ? visualSelectionRef.current : new Set(selectedResidues);
  const gradientColors = colorScales.prolint.join(', ');

  return (
    <Fade in timeout={300}>
      <Box onMouseLeave={onMouseLeave}>
        {chainIds.map((chainId) => {
          const chainResidues = residuesByChain[chainId];
          return (
            <Box key={chainId} sx={{ mb: 1.5 }}>
              <Typography variant="caption" sx={{ fontWeight: 600, mb: 0.5, display: 'block', color: colors.text.secondary }}>
                Chain {chainId} ({chainResidues.length} residues)
              </Typography>
              <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: '1.5px' }}>
                {chainResidues.map((r) => (
                  <ResidueBox
                    key={r.resid}
                    resid={r.resid}
                    resname={r.resname}
                    chainID={r.chainID || 'A'}
                    normalizedValue={normalizeValue(r.value)}
                    isSelected={currentSelection.has(r.resid)}
                    metricLabel={metricLabel}
                    rawValue={r.value}
                    onMouseDown={onMouseDown}
                    onMouseEnter={onMouseEnter}
                    onMouseUp={onMouseUp}
                  />
                ))}
              </Box>
            </Box>
          );
        })}
        <Box sx={{ mt: 1.5, display: 'flex', justifyContent: 'center', alignItems: 'center', gap: 1 }}>
          <Typography variant="caption" sx={{ fontSize: 9 }}>{minVal.toFixed(2)}</Typography>
          <Box sx={{ width: 150, height: 10, background: `linear-gradient(to right, ${gradientColors})`, borderRadius: '2px', border: `1px solid ${colors.border.default}` }} />
          <Typography variant="caption" sx={{ fontSize: 9 }}>{maxVal.toFixed(2)}</Typography>
          <Typography variant="caption" sx={{ fontSize: 9, color: colors.text.secondary, ml: 0.5 }}>{metricLabel}</Typography>
        </Box>
      </Box>
    </Fade>
  );
});

export default LogoPlotSection;
