/** Hook for logo plot visualization and residue selection with click/drag support. */

import { useState, useRef, useMemo, useCallback } from 'react';
import axios, { CancelTokenSource } from 'axios';
import type { LogoPlotData } from '../types';

interface UseLogoPlotDataProps {
  resultId: string;
  selectedMetric: string;
}

export function useLogoPlotData({ resultId, selectedMetric }: UseLogoPlotDataProps) {
  const [logoPlotData, setLogoPlotData] = useState<LogoPlotData | null>(null);
  const [logoPlotLoading, setLogoPlotLoading] = useState(false);
  const [selectedResidues, setSelectedResidues] = useState<number[]>([]);

  // Drag state refs
  const dragStartRef = useRef<number | null>(null);
  const isDraggingRef = useRef(false);
  const selectionBeforeDragRef = useRef<number[]>([]);
  const visualSelectionRef = useRef<Set<number>>(new Set());
  const [, forceUpdate] = useState(0);

  // Cancel token
  const logoPlotCancelRef = useRef<CancelTokenSource | null>(null);

  const loadLogoPlotData = useCallback(async (dbType: string, metric?: string, id?: string) => {
    const targetId = id || resultId;
    if (!targetId || !dbType) return;

    if (logoPlotCancelRef.current) {
      logoPlotCancelRef.current.cancel('New request initiated');
    }
    logoPlotCancelRef.current = axios.CancelToken.source();

    const metricToUse = metric ?? selectedMetric;

    setLogoPlotLoading(true);
    try {
      const params = { metric: metricToUse, database_type: dbType };
      const response = await axios.get(`/api/dashboard/${targetId}/logoplot`, {
        params,
        cancelToken: logoPlotCancelRef.current.token
      });
      setLogoPlotData(response.data);
    } catch (err: unknown) {
      if (!axios.isCancel(err)) {
        setLogoPlotData(null);
      }
    } finally {
      setLogoPlotLoading(false);
    }
  }, [resultId, selectedMetric]);

  // Logo plot selection handlers
  const handleResidueMouseDown = useCallback((residueId: number) => {
    dragStartRef.current = residueId;
    selectionBeforeDragRef.current = [...selectedResidues];
    visualSelectionRef.current = new Set(selectedResidues);
  }, [selectedResidues]);

  const handleResidueMouseEnter = useCallback((residueId: number) => {
    if (dragStartRef.current === null) return;
    if (residueId !== dragStartRef.current) {
      isDraggingRef.current = true;
    }

    if (!logoPlotData) return;
    const allResidueIds = logoPlotData.residues.map(r => r.resid);

    const startIdx = allResidueIds.indexOf(dragStartRef.current);
    const endIdx = allResidueIds.indexOf(residueId);
    if (startIdx === -1 || endIdx === -1) return;

    const minIdx = Math.min(startIdx, endIdx);
    const maxIdx = Math.max(startIdx, endIdx);
    const rangeResidues = allResidueIds.slice(minIdx, maxIdx + 1);

    visualSelectionRef.current = new Set([...selectionBeforeDragRef.current, ...rangeResidues]);
    forceUpdate(n => n + 1);
  }, [logoPlotData]);

  const handleResidueMouseUp = useCallback((residueId: number) => {
    if (isDraggingRef.current) {
      setSelectedResidues(Array.from(visualSelectionRef.current));
    } else if (dragStartRef.current === residueId) {
      const newSelection = selectedResidues.includes(residueId)
        ? selectedResidues.filter(r => r !== residueId)
        : [...selectedResidues, residueId];
      setSelectedResidues(newSelection);
      visualSelectionRef.current = new Set(newSelection);
    }
    isDraggingRef.current = false;
    dragStartRef.current = null;
    selectionBeforeDragRef.current = [];
  }, [selectedResidues]);

  const handleResidueMouseLeave = useCallback(() => {
    if (isDraggingRef.current && visualSelectionRef.current.size > 0) {
      setSelectedResidues(Array.from(visualSelectionRef.current));
    }
    isDraggingRef.current = false;
    dragStartRef.current = null;
    selectionBeforeDragRef.current = [];
  }, []);

  // Logo plot normalization
  const logoPlotNormalization = useMemo(() => {
    if (!logoPlotData || logoPlotData.residues.length === 0) {
      return { minVal: 0, maxVal: 1, normalizeValue: (_val: number) => 0.5 };
    }
    const values = logoPlotData.residues.map(r => r.value);
    const minVal = Math.min(...values);
    const maxVal = Math.max(...values);
    const normalizeValue = (val: number) => (maxVal === minVal) ? 0.5 : (val - minVal) / (maxVal - minVal);
    return { minVal, maxVal, normalizeValue };
  }, [logoPlotData]);

  // Residues by chain
  const residuesByChain = useMemo(() => {
    if (!logoPlotData) return {};
    return logoPlotData.residues.reduce((acc, r) => {
      const chain = r.chainID || 'A';
      if (!acc[chain]) acc[chain] = [];
      acc[chain].push(r);
      return acc;
    }, {} as Record<string, typeof logoPlotData.residues>);
  }, [logoPlotData]);

  const clearLogoPlotData = useCallback(() => {
    setLogoPlotData(null);
    setSelectedResidues([]);
  }, []);

  return {
    logoPlotData,
    logoPlotLoading,
    selectedResidues,
    setSelectedResidues,
    logoPlotNormalization,
    residuesByChain,
    isDraggingRef,
    visualSelectionRef,
    loadLogoPlotData,
    handleResidueMouseDown,
    handleResidueMouseEnter,
    handleResidueMouseUp,
    handleResidueMouseLeave,
    clearLogoPlotData,
  };
}
