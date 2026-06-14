/** Hook for 2D density map and 3D structure projection data. */

import { useState, useRef, useCallback } from 'react';
import axios, { CancelTokenSource } from 'axios';
import type { DensityMapData } from '../types';

interface UseDensityDataProps {
  resultId: string;
  startFrame: number;
  selectedMetric: string;
}

export function useDensityData({ resultId, startFrame, selectedMetric }: UseDensityDataProps) {
  const [densityMap, setDensityMap] = useState<DensityMapData | null>(null);
  const [densityMapLoading, setDensityMapLoading] = useState(false);
  const [threeDProjection, set3DProjection] = useState<string | null>(null);
  const [threeDProjectionLoading, set3DProjectionLoading] = useState(false);
  const [structureError, setStructureError] = useState<string | null>(null);

  const densityMapCancelRef = useRef<CancelTokenSource | null>(null);
  const threeDProjectionCancelRef = useRef<CancelTokenSource | null>(null);

  const loadDensityMap = useCallback(async (frameStart: number, frameEnd: number, dbType: string, id?: string) => {
    const targetId = id || resultId;
    if (!targetId || !dbType) return;

    if (densityMapCancelRef.current) {
      densityMapCancelRef.current.cancel('New request initiated');
    }
    densityMapCancelRef.current = axios.CancelToken.source();

    setDensityMapLoading(true);
    try {
      const params: Record<string, unknown> = { frame_start: frameStart, frame_end: frameEnd, bins: 150 };
      params.database_types = dbType;
      const response = await axios.get(`/api/dashboard/${targetId}/density-map`, {
        params,
        cancelToken: densityMapCancelRef.current.token,
        timeout: 30000,
      });
      setDensityMap(response.data);
    } catch (err: unknown) {
      if (!axios.isCancel(err)) {
        const axiosError = err as { message?: string };
        setStructureError(axiosError.message || 'Failed to load density map');
      }
    } finally {
      setDensityMapLoading(false);
    }
  }, [resultId]);

  const load3DProjection = useCallback(async (dbType: string, frame?: number, metric?: string, id?: string) => {
    const targetId = id || resultId;
    if (!targetId || !dbType) return;

    // Cancel any pending 3D projection request
    if (threeDProjectionCancelRef.current) {
      threeDProjectionCancelRef.current.cancel('New 3D projection request initiated');
    }
    threeDProjectionCancelRef.current = axios.CancelToken.source();

    const frameIdx = frame ?? startFrame;
    const metricToUse = metric ?? selectedMetric;

    set3DProjectionLoading(true);
    try {
      const params: Record<string, unknown> = { metric: metricToUse, frame_idx: frameIdx, database_type: dbType };
      const response = await axios.get(`/api/dashboard/${targetId}/structure`, {
        params,
        responseType: 'text',
        cancelToken: threeDProjectionCancelRef.current.token,
        timeout: 30000,
      });
      set3DProjection(response.data);
    } catch (err: unknown) {
      if (!axios.isCancel(err)) {
        set3DProjection(null);
      }
    } finally {
      set3DProjectionLoading(false);
    }
  }, [resultId, startFrame, selectedMetric]);

  const handleDownloadStructure = useCallback(async (selectedDatabaseType: string | null, endFrame: number) => {
    if (!selectedDatabaseType) return;

    setStructureError(null);
    try {
      const params: Record<string, unknown> = { metric: selectedMetric, frame_idx: startFrame, database_type: selectedDatabaseType };

      const response = await axios.get(`/api/dashboard/${resultId}/structure`, {
        params,
        responseType: 'blob'
      });

      const blob = new Blob([response.data], { type: 'chemical/x-pdb' });
      const url = window.URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = `query_structure_${selectedMetric}_${selectedDatabaseType}_frames${startFrame}-${endFrame}.pdb`;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      window.URL.revokeObjectURL(url);
    } catch (err: unknown) {
      const axiosError = err as { response?: { data?: { detail?: string } }; message?: string };
      setStructureError(axiosError.response?.data?.detail || axiosError.message || 'Failed to download structure');
    }
  }, [resultId, startFrame, selectedMetric]);

  const clearDensityData = useCallback(() => {
    setDensityMap(null);
    set3DProjection(null);
  }, []);

  return {
    densityMap,
    densityMapLoading,
    threeDProjection,
    threeDProjectionLoading,
    structureError,
    loadDensityMap,
    load3DProjection,
    handleDownloadStructure,
    clearDensityData,
  };
}
