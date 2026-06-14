/** ResidueBox component for logo plot with accessibility support. */

import { memo, useCallback } from 'react';
import { Box } from '@mui/material';
import { colors } from '../../../theme/visualizationTheme';
import { AA_MAPPING, getLogoColor } from '../utils';
import { getContrastTextColor } from '../../../utils/colorUtils';

export interface ResidueBoxProps {
  resid: number;
  resname: string;
  chainID: string;
  normalizedValue: number;
  isSelected: boolean;
  metricLabel: string;
  rawValue: number;
  onMouseDown: (resid: number) => void;
  onMouseEnter: (resid: number) => void;
  onMouseUp: (resid: number) => void;
}

/** Get WCAG-compliant text color for background. */
function getTextColor(bgColor: string, normalizedValue: number): string {
  // Try to use the color utility for hex colors
  if (bgColor.startsWith('#')) {
    return getContrastTextColor(bgColor) === 'dark' ? '#000' : '#fff';
  }

  // For rgb/rgba colors, parse and check luminance
  const rgbMatch = bgColor.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)/);
  if (rgbMatch) {
    const r = parseInt(rgbMatch[1]);
    const g = parseInt(rgbMatch[2]);
    const b = parseInt(rgbMatch[3]);
    // Simple luminance check
    const luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255;
    return luminance > 0.5 ? '#000' : '#fff';
  }

  // Fallback to original logic
  return normalizedValue > 0.5 ? '#000' : '#fff';
}

const ResidueBox = memo(function ResidueBox({
  resid,
  resname,
  chainID,
  normalizedValue,
  isSelected,
  metricLabel,
  rawValue,
  onMouseDown,
  onMouseEnter,
  onMouseUp,
}: ResidueBoxProps) {
  const oneLetter = AA_MAPPING[resname] || resname.charAt(0);
  const bgColor = getLogoColor(normalizedValue);
  const textColor = getTextColor(bgColor, normalizedValue);

  // Handle keyboard interaction for accessibility
  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        onMouseDown(resid);
        onMouseUp(resid);
      }
    },
    [resid, onMouseDown, onMouseUp]
  );

  const ariaLabel = `${resname} residue ${resid}, Chain ${chainID}, ${metricLabel}: ${rawValue.toFixed(3)}${isSelected ? ', selected' : ''}`;

  return (
    <Box
      role="button"
      tabIndex={0}
      aria-label={ariaLabel}
      aria-pressed={isSelected}
      onMouseDown={() => onMouseDown(resid)}
      onMouseEnter={() => onMouseEnter(resid)}
      onMouseUp={() => onMouseUp(resid)}
      onKeyDown={handleKeyDown}
      sx={{
        width: 25,
        height: 25,
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        justifyContent: 'center',
        backgroundColor: bgColor,
        border: `1px solid ${colors.border.default}`,
        outline: isSelected ? `2px solid ${colors.data.highlight}` : 'none',
        outlineOffset: '-1px',
        borderRadius: '2px',
        cursor: 'pointer',
        fontFamily: 'monospace',
        fontWeight: 'bold',
        fontSize: '10px',
        color: textColor,
        transition: 'transform 0.1s, box-shadow 0.1s',
        userSelect: 'none',
        '&:hover': {
          transform: 'scale(1.2)',
          zIndex: 10,
        },
        '&:focus': {
          boxShadow: `0 0 0 2px ${colors.primary[500]}`,
          zIndex: 11,
        },
        '&:focus:not(:focus-visible)': {
          boxShadow: 'none',
        },
        '&:focus-visible': {
          boxShadow: `0 0 0 2px ${colors.primary[500]}`,
          zIndex: 11,
        },
      }}
      title={`${resname} ${resid} (Chain ${chainID})\n${metricLabel}: ${rawValue.toFixed(3)}`}
    >
      <div aria-hidden="true">{oneLetter}</div>
      <div aria-hidden="true" style={{ fontSize: '7px', fontWeight: 'normal', marginTop: '-2px' }}>{resid}</div>
    </Box>
  );
});

export default ResidueBox;
