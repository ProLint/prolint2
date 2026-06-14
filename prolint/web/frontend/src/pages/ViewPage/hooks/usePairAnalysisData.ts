/** Hook for residue pair analysis: distance-time, 3D structure, atom distances, and kinetics. */

import { useState, useCallback } from 'react';
import axios from 'axios';
import type { DistanceTimeData, AtomDistancesData, KineticsData } from '../types';

interface UsePairAnalysisDataProps {
  resultId: string;
  selectedDatabaseType: string | null;
}

export function usePairAnalysisData({ resultId, selectedDatabaseType }: UsePairAnalysisDataProps) {
  const [selectedDatabaseResidue, setSelectedDatabaseResidue] = useState<number | null>(null);
  const [distanceTimeData, setDistanceTimeData] = useState<DistanceTimeData | null>(null);
  const [distanceTimeLoading, setDistanceTimeLoading] = useState(false);
  const [distanceType, setDistanceType] = useState<'com' | 'min'>('min');
  const [pairFrame, setPairFrame] = useState(0);
  const [pairStructure, setPairStructure] = useState<string | null>(null);
  const [pairStructureLoading, setPairStructureLoading] = useState(false);
  const [atomDistancesData, setAtomDistancesData] = useState<AtomDistancesData | null>(null);
  const [atomDistancesLoading, setAtomDistancesLoading] = useState(false);
  const [kineticsData, setKineticsData] = useState<KineticsData | null>(null);
  const [kineticsLoading, setKineticsLoading] = useState(false);
  const [kineticsMode, setKineticsMode] = useState<'individual' | 'accumulated'>('individual');

  const loadPairStructure = useCallback(async (queryResidue: number, dbResidue: number, frameIdx: number, id?: string) => {
    const targetId = id || resultId;
    if (!targetId || queryResidue === null || dbResidue === null) {
      setPairStructure(null);
      return;
    }

    setPairStructureLoading(true);
    try {
      const params = {
        query_residue: queryResidue,
        database_residue: dbResidue,
        frame_idx: frameIdx
      };
      const response = await axios.get(`/api/dashboard/${targetId}/structure-interaction`, {
        params,
        responseType: 'text',
      });
      setPairStructure(response.data);
    } catch {
      setPairStructure(null);
    } finally {
      setPairStructureLoading(false);
    }
  }, [resultId]);

  const loadAtomDistances = useCallback(async (queryResidue: number, dbResidue: number, frameIdx: number, id?: string) => {
    const targetId = id || resultId;
    if (!targetId || queryResidue === null || dbResidue === null) {
      setAtomDistancesData(null);
      return;
    }

    setAtomDistancesLoading(true);
    try {
      const params = {
        query_residue: queryResidue,
        database_residue: dbResidue,
        frame_idx: frameIdx
      };
      const response = await axios.get(`/api/dashboard/${targetId}/atom-distances`, { params });
      setAtomDistancesData(response.data);
    } catch {
      setAtomDistancesData(null);
    } finally {
      setAtomDistancesLoading(false);
    }
  }, [resultId]);

  const loadKineticsData = useCallback(async (
    queryResidue: number,
    dbResidue: number | null,
    dbType: string | null,
    mode: 'individual' | 'accumulated',
    id?: string
  ) => {
    const targetId = id || resultId;
    if (!targetId || queryResidue === null) {
      setKineticsData(null);
      return;
    }

    if (mode === 'individual' && dbResidue === null) {
      setKineticsData(null);
      return;
    }

    if (mode === 'accumulated' && dbType === null) {
      setKineticsData(null);
      return;
    }

    setKineticsLoading(true);
    try {
      const params: Record<string, unknown> = {
        query_residue: queryResidue,
        mode: mode
      };

      if (mode === 'individual') {
        params.database_residue = dbResidue;
      } else {
        params.database_type = dbType;
      }

      const response = await axios.get(`/api/dashboard/${targetId}/kinetics`, { params });
      setKineticsData(response.data);
    } catch {
      setKineticsData(null);
    } finally {
      setKineticsLoading(false);
    }
  }, [resultId]);

  const loadDistanceTimeData = useCallback(async (queryResidue: number, dbResidue: number, id?: string) => {
    const targetId = id || resultId;
    if (!targetId || queryResidue === null || dbResidue === null) {
      setDistanceTimeData(null);
      return;
    }

    setDistanceTimeLoading(true);
    try {
      const params = {
        query_residue: queryResidue,
        database_residue: dbResidue
      };
      const response = await axios.get(`/api/dashboard/${targetId}/distance-time`, { params });
      setDistanceTimeData(response.data);

      const frames = response.data.frames || [];
      const contactFrames = response.data.contact_frames || [];
      let initialFrameIdx = 0;
      if (contactFrames.length > 0 && frames.length > 0) {
        const contactFrameIdx = frames.indexOf(contactFrames[0]);
        if (contactFrameIdx !== -1) {
          initialFrameIdx = contactFrameIdx;
        }
      }
      setPairFrame(initialFrameIdx);
      const initialFrame = frames.length > 0 ? frames[initialFrameIdx] : 0;

      await Promise.all([
        loadPairStructure(queryResidue, dbResidue, initialFrame, targetId),
        loadAtomDistances(queryResidue, dbResidue, initialFrame, targetId),
        loadKineticsData(queryResidue, dbResidue, selectedDatabaseType, kineticsMode, targetId)
      ]);
    } catch {
      setDistanceTimeData(null);
    } finally {
      setDistanceTimeLoading(false);
    }
  }, [resultId, selectedDatabaseType, kineticsMode, loadPairStructure, loadAtomDistances, loadKineticsData]);

  // Clear all pair analysis data - defined early so it can be used by other callbacks
  const clearPairAnalysisData = useCallback(() => {
    setSelectedDatabaseResidue(null);
    setDistanceTimeData(null);
    setPairFrame(0);
    setPairStructure(null);
    setAtomDistancesData(null);
    setKineticsData(null);
  }, []);

  // Database contacts row click handler
  const handleDatabaseContactsRowClick = useCallback((rowLabel: string, selectedTimeSeriesResidue: number | null) => {
    const dbResidueId = parseInt(rowLabel, 10);
    if (isNaN(dbResidueId) || selectedTimeSeriesResidue === null) return;

    if (selectedDatabaseResidue === dbResidueId) {
      clearPairAnalysisData();
      return;
    }

    setSelectedDatabaseResidue(dbResidueId);
    loadDistanceTimeData(selectedTimeSeriesResidue, dbResidueId);
  }, [selectedDatabaseResidue, loadDistanceTimeData, clearPairAnalysisData]);

  // Pair frame select
  const handlePairFrameSelect = useCallback(async (frameIdx: number, selectedTimeSeriesResidue: number | null) => {
    if (!distanceTimeData || !selectedTimeSeriesResidue || !selectedDatabaseResidue) return;
    setPairFrame(frameIdx);
    const actualFrame = distanceTimeData.frames[frameIdx];
    await Promise.all([
      loadPairStructure(selectedTimeSeriesResidue, selectedDatabaseResidue, actualFrame),
      loadAtomDistances(selectedTimeSeriesResidue, selectedDatabaseResidue, actualFrame)
    ]);
  }, [distanceTimeData, selectedDatabaseResidue, loadPairStructure, loadAtomDistances]);

  // Kinetics mode change
  const handleKineticsModeChange = useCallback(async (newMode: 'individual' | 'accumulated', selectedTimeSeriesResidue: number | null) => {
    if (newMode === kineticsMode) return;
    setKineticsMode(newMode);

    if (selectedTimeSeriesResidue !== null) {
      await loadKineticsData(
        selectedTimeSeriesResidue,
        selectedDatabaseResidue,
        selectedDatabaseType,
        newMode
      );
    }
  }, [kineticsMode, selectedDatabaseResidue, selectedDatabaseType, loadKineticsData]);

  return {
    selectedDatabaseResidue,
    setSelectedDatabaseResidue,
    distanceTimeData,
    distanceTimeLoading,
    distanceType,
    setDistanceType,
    pairFrame,
    setPairFrame,
    pairStructure,
    pairStructureLoading,
    atomDistancesData,
    atomDistancesLoading,
    kineticsData,
    kineticsLoading,
    kineticsMode,
    loadDistanceTimeData,
    loadPairStructure,
    loadAtomDistances,
    loadKineticsData,
    handleDatabaseContactsRowClick,
    handlePairFrameSelect,
    handleKineticsModeChange,
    clearPairAnalysisData,
  };
}
