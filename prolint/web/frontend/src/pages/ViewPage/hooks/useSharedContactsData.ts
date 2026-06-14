/** Hook for shared contacts matrix and network graph visualization. */

import { useState, useRef, useMemo, useCallback } from 'react';
import axios, { CancelTokenSource } from 'axios';
import type { SharedContactsData, NetworkNode, NetworkEdge } from '../types';

interface UseSharedContactsDataProps {
  resultId: string;
  selectedResidues: number[];
}

export function useSharedContactsData({ resultId, selectedResidues }: UseSharedContactsDataProps) {
  const [sharedContactsData, setSharedContactsData] = useState<SharedContactsData | null>(null);
  const [sharedContactsLoading, setSharedContactsLoading] = useState(false);

  const sharedContactsCancelRef = useRef<CancelTokenSource | null>(null);

  const loadSharedContactsData = useCallback(async (dbType: string, id?: string) => {
    const targetId = id || resultId;
    if (!targetId || !dbType) return;

    if (sharedContactsCancelRef.current) {
      sharedContactsCancelRef.current.cancel('New request initiated');
    }
    sharedContactsCancelRef.current = axios.CancelToken.source();

    setSharedContactsLoading(true);
    try {
      const params = { database_type: dbType };
      const response = await axios.get(`/api/dashboard/${targetId}/shared-contacts`, {
        params,
        cancelToken: sharedContactsCancelRef.current.token
      });
      setSharedContactsData(response.data);
    } catch (err: unknown) {
      if (!axios.isCancel(err)) {
        setSharedContactsData(null);
      }
    } finally {
      setSharedContactsLoading(false);
    }
  }, [resultId]);

  // Network data
  const networkData = useMemo((): { nodes: NetworkNode[]; edges: NetworkEdge[] } => {
    if (!sharedContactsData || !sharedContactsData.matrix || sharedContactsData.labels.length < 2) {
      return { nodes: [], edges: [] };
    }

    if (selectedResidues.length === 0) {
      return { nodes: [], edges: [] };
    }

    const { labels, matrix } = sharedContactsData;

    const selectedSet = new Set(selectedResidues);
    const filteredIndices: number[] = [];
    const filteredLabels: number[] = [];

    labels.forEach((resid: number, idx: number) => {
      if (selectedSet.has(resid)) {
        filteredIndices.push(idx);
        filteredLabels.push(resid);
      }
    });

    const nodes: NetworkNode[] = filteredLabels.map((resid: number) => ({
      id: `res-${resid}`,
      label: `${resid}`,
      type: 'query',
      residue_id: resid,
      restype: 'query',
    }));

    const edges: NetworkEdge[] = [];
    for (let i = 0; i < filteredIndices.length; i++) {
      for (let j = i + 1; j < filteredIndices.length; j++) {
        const origI = filteredIndices[i];
        const origJ = filteredIndices[j];
        const value = matrix[origI][origJ];
        if (value > 0) {
          edges.push({
            source: `res-${filteredLabels[i]}`,
            target: `res-${filteredLabels[j]}`,
            value: value,
          });
        }
      }
    }

    return { nodes, edges };
  }, [sharedContactsData, selectedResidues]);

  const clearSharedContactsData = useCallback(() => {
    setSharedContactsData(null);
  }, []);

  return {
    sharedContactsData,
    sharedContactsLoading,
    networkData,
    loadSharedContactsData,
    clearSharedContactsData,
  };
}
