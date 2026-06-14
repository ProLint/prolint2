/** Utility functions and constants for ViewPage. */

import { useState, useEffect } from 'react';
import { colorScales } from '../../theme/visualizationTheme';

/** Debounce hook for delaying value updates. */
export function useDebounce<T>(value: T, delay: number): T {
  const [debouncedValue, setDebouncedValue] = useState<T>(value);

  useEffect(() => {
    const handler = setTimeout(() => {
      setDebouncedValue(value);
    }, delay);

    return () => {
      clearTimeout(handler);
    };
  }, [value, delay]);

  return debouncedValue;
}

export const unitLabels: Record<string, string> = {
  fs: 'fs',
  ps: 'ps',
  ns: 'ns',
  us: 'μs',
  ms: 'ms',
  s: 's',
};

export const getUnitLabel = (units: string): string => unitLabels[units] || units;

export const AA_MAPPING: Record<string, string> = {
  'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
  'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
  'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
  'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
};

export const getLogoColor = (normalizedValue: number): string => {
  const scale = colorScales.prolint;
  const idx = Math.min(Math.floor(normalizedValue * (scale.length - 1)), scale.length - 1);
  return scale[idx];
};
