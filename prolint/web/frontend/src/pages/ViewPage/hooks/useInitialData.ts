/** Hook for initial data loading, result ID parsing, and frame range management. */

import { useState, useEffect, useMemo } from 'react';
import { useSearchParams } from 'react-router-dom';
import axios from 'axios';
import { colorScales } from '../../../theme/visualizationTheme';
import type { InteractionData, CompositionItem } from '../types';

export function useInitialData() {
  const [searchParams] = useSearchParams();
  const [resultId, setResultId] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [interactionData, setInteractionData] = useState<InteractionData | null>(null);
  const [autoLoaded, setAutoLoaded] = useState(false);

  // Frame range state
  const [startFrame, setStartFrame] = useState(0);
  const [endFrame, setEndFrame] = useState(0);

  // Load initial data
  const handleLoadData = async (id?: string) => {
    const idToLoad = id || resultId;
    if (!idToLoad) {
      setError('Please enter a result ID');
      return;
    }

    if (id) {
      setResultId(id);
    }

    setLoading(true);
    setError(null);

    try {
      const response = await axios.get(`/api/dashboard/${idToLoad}/interactions`);
      setInteractionData(response.data.interactions);

      const frameRange = response.data.interactions?.frame_range;
      if (frameRange) {
        setStartFrame(frameRange.start);
        setEndFrame(frameRange.end);
      } else {
        const nFrames = response.data.interactions?.universe?.n_frames || 1;
        setStartFrame(0);
        setEndFrame(nFrames - 1);
      }
    } catch (err: unknown) {
      const axiosError = err as { response?: { data?: { detail?: string } }; message?: string };
      setError(axiosError.response?.data?.detail || axiosError.message || 'Failed to load data');
    } finally {
      setLoading(false);
    }
  };

  // Auto-load from URL
  useEffect(() => {
    const urlResultId = searchParams.get('result_id');
    if (urlResultId && !autoLoaded && !interactionData) {
      setAutoLoaded(true);
      handleLoadData(urlResultId);
    }
  }, [searchParams, autoLoaded, interactionData]);

  // Composition data
  const compositionData = useMemo((): CompositionItem[] => {
    if (!interactionData) return [];

    const resnameCounts = interactionData.composition.resname_counts;
    if (!resnameCounts || Object.keys(resnameCounts).length === 0) return [];

    const categoricalColors = colorScales.categorical;

    const typeCounts = Object.entries(resnameCounts).map(([type, count]) => ({
      type,
      count,
    }));

    const totalCount = typeCounts.reduce((sum, item) => sum + item.count, 0);
    if (totalCount === 0) return [];

    return typeCounts
      .map((item, index) => ({
        name: item.type,
        value: item.count,
        percentage: (item.count / totalCount) * 100,
        color: categoricalColors[index % categoricalColors.length],
      }))
      .sort((a, b) => b.value - a.value);
  }, [interactionData]);

  return {
    resultId,
    setResultId,
    loading,
    error,
    setError,
    interactionData,
    startFrame,
    setStartFrame,
    endFrame,
    setEndFrame,
    compositionData,
    handleLoadData,
  };
}
