/** Main composition hook combining all ViewPage domain hooks. */

import { useState, useRef, useCallback } from 'react';
import { useInitialData } from './useInitialData';
import { useDensityData } from './useDensityData';
import { useLogoPlotData } from './useLogoPlotData';
import { useSharedContactsData } from './useSharedContactsData';
import { useTimeSeriesData } from './useTimeSeriesData';
import { usePairAnalysisData } from './usePairAnalysisData';
import type { ExpandedSections } from '../types';

export function useViewPageData() {
  // Database type selection
  const [selectedDatabaseType, setSelectedDatabaseType] = useState<string | null>(null);
  const [selectedMetric, setSelectedMetric] = useState<string>('occupancy');

  // Collapsible sections
  const [expandedSections, setExpandedSections] = useState<ExpandedSections>({
    overview: true,
    residueSelection: true,
    timeSeries: true,
    pairAnalysis: true,
  });

  // Section refs
  const overviewRef = useRef<HTMLDivElement>(null);
  const residueSelectionRef = useRef<HTMLDivElement>(null);
  const timeSeriesRef = useRef<HTMLDivElement>(null);
  const pairAnalysisRef = useRef<HTMLDivElement>(null);

  // Initialize domain hooks
  const initialData = useInitialData();

  const densityData = useDensityData({
    resultId: initialData.resultId,
    startFrame: initialData.startFrame,
    selectedMetric,
  });

  const logoPlotData = useLogoPlotData({
    resultId: initialData.resultId,
    selectedMetric,
  });

  const sharedContactsData = useSharedContactsData({
    resultId: initialData.resultId,
    selectedResidues: logoPlotData.selectedResidues,
  });

  const timeSeriesData = useTimeSeriesData({
    resultId: initialData.resultId,
    selectedResidues: logoPlotData.selectedResidues,
    selectedDatabaseType,
  });

  const pairAnalysisData = usePairAnalysisData({
    resultId: initialData.resultId,
    selectedDatabaseType,
  });

  // Toggle section
  const toggleSection = useCallback((section: keyof ExpandedSections) => {
    setExpandedSections(prev => ({ ...prev, [section]: !prev[section] }));
  }, []);

  // Navigate to section
  const navigateToSection = useCallback((sectionKey: keyof ExpandedSections) => {
    const refs: Record<string, React.RefObject<HTMLDivElement>> = {
      overview: overviewRef,
      residueSelection: residueSelectionRef,
      timeSeries: timeSeriesRef,
      pairAnalysis: pairAnalysisRef,
    };
    const ref = refs[sectionKey];

    if (!expandedSections[sectionKey]) {
      setExpandedSections(prev => ({ ...prev, [sectionKey]: true }));
    }
    setTimeout(() => {
      ref?.current?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }, 100);
  }, [expandedSections]);

  // Handle database type click - fire requests without blocking
  const handleDatabaseTypeClick = useCallback((typeName: string) => {
    if (selectedDatabaseType === typeName) {
      // Deselect
      setSelectedDatabaseType(null);
      densityData.clearDensityData();
      logoPlotData.clearLogoPlotData();
      sharedContactsData.clearSharedContactsData();
      timeSeriesData.clearTimeSeriesData();
      pairAnalysisData.clearPairAnalysisData();
      return;
    }

    // Select new type - state updates immediately
    setSelectedDatabaseType(typeName);
    logoPlotData.setSelectedResidues([]);
    timeSeriesData.clearTimeSeriesData();
    pairAnalysisData.clearPairAnalysisData();

    // Fire all requests in parallel WITHOUT awaiting - each updates UI independently
    densityData.loadDensityMap(initialData.startFrame, initialData.endFrame, typeName);
    densityData.load3DProjection(typeName, initialData.startFrame, selectedMetric);
    logoPlotData.loadLogoPlotData(typeName, selectedMetric);
    sharedContactsData.loadSharedContactsData(typeName);
  }, [
    selectedDatabaseType, initialData.startFrame, initialData.endFrame, selectedMetric,
    densityData.clearDensityData, densityData.loadDensityMap, densityData.load3DProjection,
    logoPlotData.clearLogoPlotData, logoPlotData.setSelectedResidues, logoPlotData.loadLogoPlotData,
    sharedContactsData.clearSharedContactsData, sharedContactsData.loadSharedContactsData,
    timeSeriesData.clearTimeSeriesData, pairAnalysisData.clearPairAnalysisData
  ]);

  // Handle frame range change
  const handleFrameRangeChange = useCallback(async () => {
    if (!selectedDatabaseType) return;
    await Promise.all([
      densityData.loadDensityMap(initialData.startFrame, initialData.endFrame, selectedDatabaseType),
      densityData.load3DProjection(selectedDatabaseType, initialData.startFrame, selectedMetric)
    ]);
  }, [selectedDatabaseType, initialData.startFrame, initialData.endFrame, selectedMetric, densityData.loadDensityMap, densityData.load3DProjection]);

  // Handle metric change
  const handleMetricChange = useCallback(async (newMetric: string) => {
    setSelectedMetric(newMetric);
    if (!selectedDatabaseType) return;
    await Promise.all([
      densityData.load3DProjection(selectedDatabaseType, initialData.startFrame, newMetric),
      logoPlotData.loadLogoPlotData(selectedDatabaseType, newMetric)
    ]);
  }, [selectedDatabaseType, initialData.startFrame, densityData.load3DProjection, logoPlotData.loadLogoPlotData]);

  // Handle download structure
  const handleDownloadStructure = useCallback(async () => {
    await densityData.handleDownloadStructure(selectedDatabaseType, initialData.endFrame);
  }, [densityData.handleDownloadStructure, selectedDatabaseType, initialData.endFrame]);

  // Handle time series row click
  const handleTimeSeriesRowClick = useCallback((rowLabel: string) => {
    timeSeriesData.handleTimeSeriesRowClick(rowLabel, () => {
      pairAnalysisData.clearPairAnalysisData();
    });
  }, [timeSeriesData.handleTimeSeriesRowClick, pairAnalysisData.clearPairAnalysisData]);

  // Handle database contacts row click
  const handleDatabaseContactsRowClick = useCallback((rowLabel: string) => {
    pairAnalysisData.handleDatabaseContactsRowClick(rowLabel, timeSeriesData.selectedTimeSeriesResidue);
  }, [pairAnalysisData.handleDatabaseContactsRowClick, timeSeriesData.selectedTimeSeriesResidue]);

  // Handle pair frame select
  const handlePairFrameSelect = useCallback(async (frameIdx: number) => {
    await pairAnalysisData.handlePairFrameSelect(frameIdx, timeSeriesData.selectedTimeSeriesResidue);
  }, [pairAnalysisData.handlePairFrameSelect, timeSeriesData.selectedTimeSeriesResidue]);

  // Handle kinetics mode change
  const handleKineticsModeChange = useCallback(async (newMode: 'individual' | 'accumulated') => {
    await pairAnalysisData.handleKineticsModeChange(newMode, timeSeriesData.selectedTimeSeriesResidue);
  }, [pairAnalysisData.handleKineticsModeChange, timeSeriesData.selectedTimeSeriesResidue]);

  // Clear time series selection
  const clearTimeSeriesSelection = useCallback(() => {
    timeSeriesData.clearTimeSeriesSelection();
    pairAnalysisData.clearPairAnalysisData();
  }, [timeSeriesData.clearTimeSeriesSelection, pairAnalysisData.clearPairAnalysisData]);

  return {
    // Initial data
    resultId: initialData.resultId,
    setResultId: initialData.setResultId,
    loading: initialData.loading,
    error: initialData.error,
    setError: initialData.setError,
    interactionData: initialData.interactionData,
    startFrame: initialData.startFrame,
    setStartFrame: initialData.setStartFrame,
    endFrame: initialData.endFrame,
    setEndFrame: initialData.setEndFrame,
    compositionData: initialData.compositionData,
    handleLoadData: initialData.handleLoadData,

    // Density data
    densityMap: densityData.densityMap,
    densityMapLoading: densityData.densityMapLoading,
    threeDProjection: densityData.threeDProjection,
    threeDProjectionLoading: densityData.threeDProjectionLoading,
    structureError: densityData.structureError,

    // Logo plot data
    logoPlotData: logoPlotData.logoPlotData,
    logoPlotLoading: logoPlotData.logoPlotLoading,
    selectedResidues: logoPlotData.selectedResidues,
    setSelectedResidues: logoPlotData.setSelectedResidues,
    logoPlotNormalization: logoPlotData.logoPlotNormalization,
    residuesByChain: logoPlotData.residuesByChain,
    isDraggingRef: logoPlotData.isDraggingRef,
    visualSelectionRef: logoPlotData.visualSelectionRef,
    handleResidueMouseDown: logoPlotData.handleResidueMouseDown,
    handleResidueMouseEnter: logoPlotData.handleResidueMouseEnter,
    handleResidueMouseUp: logoPlotData.handleResidueMouseUp,
    handleResidueMouseLeave: logoPlotData.handleResidueMouseLeave,

    // Shared contacts data
    sharedContactsData: sharedContactsData.sharedContactsData,
    sharedContactsLoading: sharedContactsData.sharedContactsLoading,
    networkData: sharedContactsData.networkData,

    // Time series data
    timeSeriesData: timeSeriesData.timeSeriesData,
    timeSeriesLoading: timeSeriesData.timeSeriesLoading,
    selectedTimeSeriesResidue: timeSeriesData.selectedTimeSeriesResidue,
    residueTimeSeriesData: timeSeriesData.residueTimeSeriesData,
    residueTimeSeriesLoading: timeSeriesData.residueTimeSeriesLoading,

    // Pair analysis data
    selectedDatabaseResidue: pairAnalysisData.selectedDatabaseResidue,
    distanceTimeData: pairAnalysisData.distanceTimeData,
    distanceTimeLoading: pairAnalysisData.distanceTimeLoading,
    distanceType: pairAnalysisData.distanceType,
    setDistanceType: pairAnalysisData.setDistanceType,
    pairFrame: pairAnalysisData.pairFrame,
    pairStructure: pairAnalysisData.pairStructure,
    pairStructureLoading: pairAnalysisData.pairStructureLoading,
    atomDistancesData: pairAnalysisData.atomDistancesData,
    atomDistancesLoading: pairAnalysisData.atomDistancesLoading,
    kineticsData: pairAnalysisData.kineticsData,
    kineticsLoading: pairAnalysisData.kineticsLoading,
    kineticsMode: pairAnalysisData.kineticsMode,

    // Selection state
    selectedDatabaseType,
    selectedMetric,
    expandedSections,

    // Refs
    overviewRef,
    residueSelectionRef,
    timeSeriesRef,
    pairAnalysisRef,

    // Handlers
    handleDownloadStructure,
    handleFrameRangeChange,
    handleMetricChange,
    handleDatabaseTypeClick,
    handleTimeSeriesRowClick,
    handleDatabaseContactsRowClick,
    handlePairFrameSelect,
    handleKineticsModeChange,
    toggleSection,
    navigateToSection,
    clearTimeSeriesSelection,
  };
}
