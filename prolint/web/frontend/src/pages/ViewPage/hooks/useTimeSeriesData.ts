/** Hook for time series heatmap data and per-residue analysis. */

import { useState, useRef, useEffect, useCallback } from 'react';
import axios, { CancelTokenSource } from 'axios';
import { useDebounce } from '../utils';
import { TIMEOUTS } from '../../../config/constants';
import type { TimeSeriesData, ResidueTimeSeriesData } from '../types';

interface UseTimeSeriesDataProps {
  resultId: string;
  selectedResidues: number[];
  selectedDatabaseType: string | null;
}

export function useTimeSeriesData({ resultId, selectedResidues, selectedDatabaseType }: UseTimeSeriesDataProps) {
  const [timeSeriesData, setTimeSeriesData] = useState<TimeSeriesData | null>(null);
  const [timeSeriesLoading, setTimeSeriesLoading] = useState(false);
  const [selectedTimeSeriesResidue, setSelectedTimeSeriesResidue] = useState<number | null>(null);
  const [residueTimeSeriesData, setResidueTimeSeriesData] = useState<ResidueTimeSeriesData | null>(null);
  const [residueTimeSeriesLoading, setResidueTimeSeriesLoading] = useState(false);

  const timeSeriesCancelRef = useRef<CancelTokenSource | null>(null);
  const residueTimeSeriesCancelRef = useRef<CancelTokenSource | null>(null);

  // Debounced selected residues for time series loading
  const debouncedSelectedResidues = useDebounce(selectedResidues, TIMEOUTS.RESIDUE_SELECTION_DEBOUNCE_MS);

  const loadTimeSeriesData = useCallback(async (queryResidues: number[], dbType: string, id?: string) => {
    const targetId = id || resultId;
    if (!targetId || !dbType || queryResidues.length === 0) {
      setTimeSeriesData(null);
      return;
    }

    if (timeSeriesCancelRef.current) {
      timeSeriesCancelRef.current.cancel('New request initiated');
    }
    timeSeriesCancelRef.current = axios.CancelToken.source();

    setTimeSeriesLoading(true);
    try {
      const params = {
        database_type: dbType,
        query_residues: queryResidues.join(',')
      };
      const response = await axios.get(`/api/dashboard/${targetId}/timeseries`, {
        params,
        cancelToken: timeSeriesCancelRef.current.token,
        timeout: 30000,
      });
      setTimeSeriesData(response.data);
    } catch (err: unknown) {
      if (!axios.isCancel(err)) {
        setTimeSeriesData(null);
      }
    } finally {
      setTimeSeriesLoading(false);
    }
  }, [resultId]);

  const loadResidueTimeSeriesData = useCallback(async (queryResidue: number, dbType: string, id?: string) => {
    const targetId = id || resultId;
    if (!targetId || !dbType || queryResidue === null) {
      setResidueTimeSeriesData(null);
      return;
    }

    // Cancel any pending residue time series request
    if (residueTimeSeriesCancelRef.current) {
      residueTimeSeriesCancelRef.current.cancel('New residue time series request initiated');
    }
    residueTimeSeriesCancelRef.current = axios.CancelToken.source();

    setResidueTimeSeriesLoading(true);
    try {
      const params = {
        query_residue: queryResidue,
        database_type: dbType
      };
      const response = await axios.get(`/api/dashboard/${targetId}/residue-timeseries`, {
        params,
        cancelToken: residueTimeSeriesCancelRef.current.token,
        timeout: 30000,
      });
      setResidueTimeSeriesData(response.data);
    } catch (err: unknown) {
      if (!axios.isCancel(err)) {
        setResidueTimeSeriesData(null);
      }
    } finally {
      setResidueTimeSeriesLoading(false);
    }
  }, [resultId]);

  // Load time series when selection changes
  useEffect(() => {
    if (debouncedSelectedResidues.length > 0 && selectedDatabaseType) {
      loadTimeSeriesData(debouncedSelectedResidues, selectedDatabaseType);
    } else {
      setTimeSeriesData(null);
    }
  }, [debouncedSelectedResidues, selectedDatabaseType, loadTimeSeriesData]);

  // Time series row click handler
  const handleTimeSeriesRowClick = useCallback((rowLabel: string, onResidueDeselected: () => void) => {
    const residueId = parseInt(rowLabel, 10);
    if (isNaN(residueId)) return;

    if (selectedTimeSeriesResidue === residueId) {
      setSelectedTimeSeriesResidue(null);
      setResidueTimeSeriesData(null);
      onResidueDeselected();
      return;
    }

    setSelectedTimeSeriesResidue(residueId);
    if (selectedDatabaseType) {
      loadResidueTimeSeriesData(residueId, selectedDatabaseType);
    }
  }, [selectedTimeSeriesResidue, selectedDatabaseType, loadResidueTimeSeriesData]);

  const clearTimeSeriesData = useCallback(() => {
    setTimeSeriesData(null);
    setSelectedTimeSeriesResidue(null);
    setResidueTimeSeriesData(null);
  }, []);

  const clearTimeSeriesSelection = useCallback(() => {
    setSelectedTimeSeriesResidue(null);
    setResidueTimeSeriesData(null);
  }, []);

  return {
    timeSeriesData,
    timeSeriesLoading,
    selectedTimeSeriesResidue,
    setSelectedTimeSeriesResidue,
    residueTimeSeriesData,
    residueTimeSeriesLoading,
    loadTimeSeriesData,
    loadResidueTimeSeriesData,
    handleTimeSeriesRowClick,
    clearTimeSeriesData,
    clearTimeSeriesSelection,
  };
}
