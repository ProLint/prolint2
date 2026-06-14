/** Pair analysis: distance plot, spatial positions, atom distances, and kinetics. */

import { memo, useMemo, useCallback } from 'react';
import { Box, Typography, Stack, Chip, Skeleton } from '@mui/material';
import ReactECharts from 'echarts-for-react';
import { colors } from '../../../theme/visualizationTheme';
import { ChartSkeleton, HeatmapSkeleton } from './Skeletons';
import { getUnitLabel } from '../utils';
import type { DistanceTimeData, DensityMapData, AtomDistancesData, KineticsData, InteractionData } from '../types';
import type { EChartsInstance, EChartsTooltipParam, ZRenderEventParams, HeatmapDataPoint } from '../../../types/echarts';

// Distance vs Time Plot
interface DistanceTimePlotProps {
  data: DistanceTimeData | null;
  loading: boolean;
  distanceType: 'com' | 'min';
  pairFrame: number;
  interactionData: InteractionData | null;
  onFrameSelect: (frameIdx: number) => void;
}

export const DistanceTimePlot = memo(function DistanceTimePlot({
  data,
  loading,
  distanceType,
  pairFrame,
  interactionData,
  onFrameSelect,
}: DistanceTimePlotProps) {
  const option = useMemo(() => {
    if (!data) return null;

    const { frames, distances, min_distances, contact_frames } = data;
    const activeDistances = distanceType === 'min' ? min_distances : distances;
    const distLabel = distanceType === 'min' ? 'Min Distance' : 'COM Distance';
    const maxDistance = Math.max(...activeDistances);

    const normalizeBy = interactionData?.params?.normalize_by || 'counts';
    const units = interactionData?.params?.units || 'ns';
    const normFactor = interactionData?.params?.norm_factor || 1;
    const showAsTime = normalizeBy === 'actual_time';
    const unitLabel = getUnitLabel(units);

    const xAxisData = showAsTime ? frames.map((f: number) => (f * normFactor).toFixed(2)) : frames.map(String);
    const contactData = contact_frames.map((frame) => {
      const idx = frames.indexOf(frame);
      const xVal = showAsTime ? (frame * normFactor).toFixed(2) : String(frame);
      return idx >= 0 ? [xVal, activeDistances[idx]] : null;
    }).filter((d): d is [string, number] => d !== null);

    const selectedFrameMarker = pairFrame >= 0 && pairFrame < frames.length
      ? [[showAsTime ? (frames[pairFrame] * normFactor).toFixed(2) : String(frames[pairFrame]), activeDistances[pairFrame]]]
      : [];

    return {
      tooltip: {
        trigger: 'axis',
        formatter: (params: EChartsTooltipParam | EChartsTooltipParam[]) => {
          const paramsArray = Array.isArray(params) ? params : [params];
          if (paramsArray.length === 0) return '';
          const xVal = paramsArray[0].axisValue;
          const dist = paramsArray[0].data;
          const frameNum = showAsTime ? Math.round(parseFloat(String(xVal)) / normFactor) : Number(xVal);
          const isContact = contact_frames.includes(frameNum);
          const timeLabel = showAsTime ? `Time: ${xVal} ${unitLabel}` : `Frame: ${xVal}`;
          const distValue = typeof dist === 'number' ? dist : Array.isArray(dist) ? (dist as number[])[1] : 0;
          return `${timeLabel}<br/>${distLabel}: ${distValue.toFixed(2)} Å${isContact ? `<br/><span style="color:${colors.success.main}">● Contact</span>` : ''}<br/><span style="color:${colors.text.tertiary}">Click to view 3D</span>`;
        },
      },
      grid: { top: 30, right: 20, bottom: 50, left: 60 },
      xAxis: {
        type: 'category',
        data: xAxisData,
        name: showAsTime ? `Time (${unitLabel})` : 'Frame',
        nameLocation: 'middle',
        nameGap: 25,
        nameTextStyle: { fontSize: 11, fontWeight: 500 },
        axisLabel: { interval: Math.max(0, Math.floor(frames.length / 10) - 1), fontSize: 9 },
      },
      yAxis: {
        type: 'value',
        name: 'Distance (Å)',
        nameLocation: 'middle',
        nameGap: 35,
        nameTextStyle: { fontSize: 11, fontWeight: 500 },
        min: 0,
        max: Math.ceil(maxDistance * 1.1),
        axisLabel: { fontSize: 9, formatter: (v: number) => v.toFixed(1) },
        splitLine: { lineStyle: { color: colors.border.light } },
      },
      series: [
        { name: distLabel, type: 'line', data: activeDistances, smooth: false, lineStyle: { color: colors.data.query, width: 2 }, areaStyle: { color: `${colors.data.query}1a` }, symbol: 'none' },
        { name: 'Contact', type: 'scatter', data: contactData, symbolSize: 5, itemStyle: { color: colors.success.main }, z: 10 },
        {
          name: 'Selected Frame', type: 'scatter', data: selectedFrameMarker, symbolSize: 14, symbol: 'diamond',
          itemStyle: { color: colors.accent[500], borderColor: colors.background.default, borderWidth: 2 }, z: 20,
          markLine: pairFrame >= 0 && pairFrame < frames.length ? {
            silent: true, symbol: 'none', lineStyle: { color: colors.accent[500], width: 2, type: 'solid' },
            data: [{ xAxis: showAsTime ? (frames[pairFrame] * normFactor).toFixed(2) : String(frames[pairFrame]) }], label: { show: false }
          } : undefined
        },
      ],
      _meta: { frames, showAsTime, normFactor }
    };
  }, [data, distanceType, pairFrame, interactionData]);

  const onChartReady = useCallback((chart: EChartsInstance) => {
    if (!chart || !option?._meta) return;
    const { frames } = option._meta;
    chart.getZr().off('click');
    chart.getZr().on('click', (params: ZRenderEventParams) => {
      const pointInPixel = [params.offsetX, params.offsetY];
      if (chart.containPixel('grid', pointInPixel)) {
        const pointInGrid = chart.convertFromPixel({ seriesIndex: 0 }, pointInPixel);
        if (pointInGrid) {
          const frameIdx = Math.round(pointInGrid[0]);
          if (frameIdx >= 0 && frameIdx < frames.length) {
            onFrameSelect(frameIdx);
          }
        }
      }
    });
  }, [option, onFrameSelect]);

  if (loading || !data) return <ChartSkeleton height={300} />;

  return <ReactECharts option={option} style={{ height: 250, width: '100%', cursor: 'crosshair' }} opts={{ renderer: 'canvas' }} notMerge={true} onChartReady={onChartReady} />;
});

// Spatial Positions Plot
interface SpatialPositionsPlotProps {
  distanceTimeData: DistanceTimeData | null;
  densityMap: DensityMapData | null;
  pairFrame: number;
  selectedTimeSeriesResidue: number | null;
  selectedDatabaseResidue: number | null;
  selectedDatabaseType: string | null;
}

export const SpatialPositionsPlot = memo(function SpatialPositionsPlot({
  distanceTimeData,
  densityMap,
  pairFrame,
  selectedTimeSeriesResidue,
  selectedDatabaseResidue,
  selectedDatabaseType,
}: SpatialPositionsPlotProps) {
  const option = useMemo(() => {
    if (!distanceTimeData || !distanceTimeData.positions) return null;

    const { positions } = distanceTimeData;
    const queryPos = positions.query[pairFrame];
    const dbPos = positions.database[pairFrame];

    if (!queryPos || !dbPos) return null;

    let minX: number, maxX: number, minY: number, maxY: number;
    if (densityMap?.x_edges && densityMap?.y_edges) {
      minX = densityMap.x_edges[0];
      maxX = densityMap.x_edges[densityMap.x_edges.length - 1];
      minY = densityMap.y_edges[0];
      maxY = densityMap.y_edges[densityMap.y_edges.length - 1];
    } else {
      const allX = [...positions.query.map(p => p.x), ...positions.database.map(p => p.x)];
      const allY = [...positions.query.map(p => p.y), ...positions.database.map(p => p.y)];
      minX = Math.min(...allX) - 10;
      maxX = Math.max(...allX) + 10;
      minY = Math.min(...allY) - 10;
      maxY = Math.max(...allY) + 10;
    }

    const densityPoints: [number, number, number][] = [];
    let maxDensity = 0;

    if (densityMap?.query_density && densityMap?.x_edges && densityMap?.y_edges) {
      densityMap.query_density.forEach((row: number[], yIdx: number) => {
        const yCenter = (densityMap.y_edges[yIdx] + densityMap.y_edges[yIdx + 1]) / 2;
        row.forEach((value: number, xIdx: number) => {
          if (value > 0) {
            const xCenter = (densityMap.x_edges[xIdx] + densityMap.x_edges[xIdx + 1]) / 2;
            densityPoints.push([xCenter, yCenter, value]);
            if (value > maxDensity) maxDensity = value;
          }
        });
      });
    }

    return {
      tooltip: { trigger: 'item', formatter: (params: EChartsTooltipParam) => {
        const d = params.data as number[];
        if (params.seriesName === 'Query Density') return `Density: ${d[2].toExponential(1)}<br/>X: ${d[0].toFixed(1)} Å<br/>Y: ${d[1].toFixed(1)} Å`;
        return `${params.seriesName}<br/>X: ${d[0].toFixed(1)} Å<br/>Y: ${d[1].toFixed(1)} Å`;
      }},
      grid: { left: 50, right: 20, top: 20, bottom: 50, containLabel: false },
      xAxis: { type: 'value', name: 'X (Å)', nameLocation: 'middle', nameGap: 25, nameTextStyle: { fontSize: 11, fontWeight: 500 }, min: minX, max: maxX, axisLabel: { fontSize: 9, formatter: (v: number) => Number(v).toFixed(0) }, splitLine: { lineStyle: { color: `${colors.text.inverse}1a` } } },
      yAxis: { type: 'value', name: 'Y (Å)', nameLocation: 'middle', nameGap: 35, nameTextStyle: { fontSize: 11, fontWeight: 500 }, min: minY, max: maxY, axisLabel: { fontSize: 9, formatter: (v: number) => Number(v).toFixed(0) }, splitLine: { lineStyle: { color: `${colors.text.inverse}1a` } } },
      visualMap: densityPoints.length > 0 ? { show: false, min: 0, max: maxDensity, dimension: 2, seriesIndex: 0, inRange: { color: [`${colors.neutral[500]}33`, `${colors.neutral[400]}80`, `${colors.neutral[300]}cc`] } } : undefined,
      backgroundColor: colors.neutral[900],
      series: [
        ...(densityPoints.length > 0 ? [{ type: 'scatter', name: 'Query Density', data: densityPoints, symbolSize: 8, itemStyle: { color: colors.neutral[400] }, z: 1 }] : []),
        { type: 'line', name: 'Distance', data: [[queryPos.x, queryPos.y], [dbPos.x, dbPos.y]], lineStyle: { color: `${colors.text.inverse}99`, width: 1, type: 'dashed' }, symbol: 'none', z: 5 },
        { type: 'scatter', data: [[queryPos.x, queryPos.y]], symbolSize: 24, itemStyle: { color: colors.data.query, borderColor: colors.text.inverse, borderWidth: 3 }, name: `Query ${selectedTimeSeriesResidue}`, z: 10 },
        { type: 'scatter', data: [[dbPos.x, dbPos.y]], symbolSize: 24, itemStyle: { color: colors.data.database, borderColor: colors.text.inverse, borderWidth: 3 }, name: `${selectedDatabaseType} ${selectedDatabaseResidue}`, z: 10 },
      ],
      _meta: { densityPoints, selectedTimeSeriesResidue, selectedDatabaseType, selectedDatabaseResidue }
    };
  }, [distanceTimeData, densityMap, pairFrame, selectedTimeSeriesResidue, selectedDatabaseResidue, selectedDatabaseType]);

  if (!option) return <ChartSkeleton height={250} />;

  const { densityPoints } = option._meta || {};

  return (
    <Box sx={{ width: '100%' }}>
      <Box sx={{ width: '100%', aspectRatio: { xs: '1 / 1.1', sm: '1 / 1' }, minHeight: { xs: 280, sm: 350 } }}>
        <ReactECharts option={option} style={{ width: '100%', height: '100%' }} opts={{ renderer: 'canvas' }} notMerge={true} />
      </Box>
      <Stack direction="row" spacing={{ xs: 1, sm: 2 }} justifyContent="center" sx={{ mt: 1, mb: 1 }} flexWrap="wrap" useFlexGap>
        {densityPoints && densityPoints.length > 0 && (
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
            <Box sx={{ width: 12, height: 12, backgroundColor: colors.neutral[400], borderRadius: '50%' }} />
            <Typography variant="caption">Query (all)</Typography>
          </Box>
        )}
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
          <Box sx={{ width: 12, height: 12, backgroundColor: colors.data.query, borderRadius: '50%' }} />
          <Typography variant="caption">Query {selectedTimeSeriesResidue}</Typography>
        </Box>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
          <Box sx={{ width: 12, height: 12, backgroundColor: colors.data.database, borderRadius: '50%' }} />
          <Typography variant="caption">{selectedDatabaseType} {selectedDatabaseResidue}</Typography>
        </Box>
      </Stack>
    </Box>
  );
});

// Atom Distance Heatmap
interface AtomDistanceHeatmapProps {
  data: AtomDistancesData | null;
  loading: boolean;
  selectedTimeSeriesResidue: number | null;
  selectedDatabaseResidue: number | null;
}

export const AtomDistanceHeatmap = memo(function AtomDistanceHeatmap({
  data,
  loading,
  selectedTimeSeriesResidue,
  selectedDatabaseResidue,
}: AtomDistanceHeatmapProps) {
  const option = useMemo(() => {
    if (!data || !data.distance_matrix) return null;

    const { query_atoms, database_atoms, distance_matrix, min_distance, max_distance } = data;
    const echartsData: [number, number, number][] = [];
    distance_matrix.forEach((row, yIdx) => {
      row.forEach((value, xIdx) => {
        echartsData.push([xIdx, yIdx, value]);
      });
    });

    return {
      tooltip: { trigger: 'item', formatter: (params: EChartsTooltipParam) => {
        const [xIdx, yIdx, dist] = params.data as HeatmapDataPoint;
        return `Query: ${query_atoms[yIdx]}<br/>Database: ${database_atoms[xIdx]}<br/>Distance: ${dist.toFixed(2)} Å`;
      }},
      grid: { left: 50, right: 10, top: 10, bottom: 60, containLabel: false },
      xAxis: { type: 'category', data: database_atoms, name: `DB ${selectedDatabaseResidue}`, nameLocation: 'middle', nameGap: 25, nameTextStyle: { fontSize: 11, fontWeight: 500 }, axisLabel: { fontSize: 9, rotate: 90, interval: 0 }, splitArea: { show: false } },
      yAxis: { type: 'category', data: query_atoms, name: `Query ${selectedTimeSeriesResidue}`, nameLocation: 'middle', nameGap: 35, nameTextStyle: { fontSize: 11, fontWeight: 500 }, axisLabel: { fontSize: 9, interval: 0 }, splitArea: { show: false }, inverse: true },
      visualMap: { min: min_distance, max: max_distance, calculable: true, orient: 'horizontal', left: 'center', bottom: -15, itemHeight: 250, itemWidth: 8, textStyle: { fontSize: 9 }, formatter: (v: number) => v.toFixed(1), inRange: { color: [colors.primary[700], colors.accent[400], colors.error.dark] } },
      series: [{ type: 'heatmap', data: echartsData, animation: false }],
    };
  }, [data, selectedTimeSeriesResidue, selectedDatabaseResidue]);

  if (loading) return <HeatmapSkeleton height={250} rows={10} />;

  if (!option) {
    return (
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: { xs: 200, sm: 250 } }}>
        <Typography color="text.secondary">No atom distance data available</Typography>
      </Box>
    );
  }

  return <ReactECharts option={option} style={{ height: '100%', width: '100%' }} opts={{ renderer: 'canvas' }} notMerge={true} />;
});

// Kinetics Metrics
interface KineticsMetricsProps {
  data: KineticsData | null;
  loading: boolean;
  interactionData: InteractionData | null;
}

export const KineticsMetrics = memo(function KineticsMetrics({
  data,
  loading,
  interactionData,
}: KineticsMetricsProps) {
  if (loading) {
    return (
      <Box sx={{ display: 'flex', gap: 0.5, justifyContent: 'center', flexWrap: 'wrap' }}>
        {[1, 2, 3, 4, 5].map((i) => <Skeleton key={i} variant="rounded" width={80} height={22} />)}
      </Box>
    );
  }

  if (!data) {
    return (
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', py: 1 }}>
        <Typography variant="caption" color="text.secondary">No kinetics data available</Typography>
      </Box>
    );
  }

  const k = data.kinetics;
  const units = interactionData?.params?.units || 'ns';
  const normalizeBy = interactionData?.params?.normalize_by || 'actual_time';

  const formatResidenceTime = (frames: number): string => {
    if (normalizeBy === 'counts') return `${frames.toFixed(1)} fr`;
    const unit = getUnitLabel(units);
    if (frames < 0.001 && frames !== 0) return frames.toExponential(2) + ' ' + unit;
    return `${frames.toFixed(2)} ${unit}`;
  };

  const formatRate = (rate: number): string => {
    if (normalizeBy === 'counts') return `${rate.toFixed(3)} fr⁻¹`;
    const unit = getUnitLabel(units);
    return `${rate.toFixed(3)} ${unit}⁻¹`;
  };

  const metrics = [
    { label: 'Occupancy', value: `${(k.occupancy * 100).toFixed(1)}%` },
    { label: 'Mean Residence', value: formatResidenceTime(k.mean_residence_time) },
    { label: 'Events', value: `${k.n_events}` },
    { label: 'k_off', value: formatRate(k.koff) },
    { label: 'K_d', value: k.kd !== null ? k.kd.toFixed(3) : '∞' },
  ];

  return (
    <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.5, justifyContent: 'center' }}>
      {metrics.map((metric, idx) => (
        <Chip key={idx} label={<><span style={{ opacity: 0.7 }}>{metric.label}:</span> <strong>{metric.value}</strong></>} size="small" variant="outlined" sx={{ fontSize: '0.65rem', height: 22 }} />
      ))}
    </Box>
  );
});

// Survival Curve
interface SurvivalCurveProps {
  data: KineticsData | null;
  interactionData: InteractionData | null;
}

export const SurvivalCurve = memo(function SurvivalCurve({
  data,
  interactionData,
}: SurvivalCurveProps) {
  const option = useMemo(() => {
    if (!data || !data.survival_curve) return null;

    const { lag_times, survival_probability, mono_fit, bi_fit, selected_model } = data.survival_curve;
    const units = interactionData?.params?.units || 'ns';
    const normalizeBy = interactionData?.params?.normalize_by || 'actual_time';
    const lagAxisLabel = normalizeBy === 'counts' ? 'Lag (fr)' : `Lag (${getUnitLabel(units)})`;

    if (lag_times.length === 0) return null;

    interface SeriesConfig {
      name: string;
      type: string;
      data: number[];
      lineStyle: { color: string; width: number };
      symbol?: string;
      symbolSize?: number;
      itemStyle?: { color: string };
    }

    const series: SeriesConfig[] = [
      { name: 'Empirical', type: 'line', data: survival_probability, lineStyle: { color: colors.data.database, width: 2 }, symbol: 'circle', symbolSize: 4, itemStyle: { color: colors.data.database } },
    ];

    if (selected_model === 'biexponential' && bi_fit) {
      series.push({ name: `Fit (k₁=${bi_fit.k_fast.toFixed(3)}, k₂=${bi_fit.k_slow.toFixed(3)})`, type: 'line', data: bi_fit.fitted_curve, lineStyle: { color: colors.data.highlight, width: 2 }, symbol: 'none' });
    } else if (selected_model === 'monoexponential' && mono_fit) {
      series.push({ name: `Fit (k=${mono_fit.k_off.toFixed(4)})`, type: 'line', data: mono_fit.fitted_curve, lineStyle: { color: colors.data.highlight, width: 2 }, symbol: 'none' });
    }

    return {
      tooltip: {
        trigger: 'axis',
        formatter: (params: EChartsTooltipParam | EChartsTooltipParam[]) => {
          const paramsArray = Array.isArray(params) ? params : [params];
          if (paramsArray.length === 0) return '';
          const lag = paramsArray[0].axisValue;
          const lagUnit = normalizeBy === 'counts' ? 'fr' : getUnitLabel(units);
          let result = `Lag: ${lag}${lagUnit ? ' ' + lagUnit : ''}`;
          paramsArray.forEach((p: EChartsTooltipParam) => {
            const val = typeof p.data === 'number' ? p.data : (p.data as number[])[1];
            result += `<br/>${p.seriesName}: ${(val * 100).toFixed(1)}%`;
          });
          return result;
        },
      },
      legend: { show: true, orient: 'horizontal', top: 0, left: 'center', textStyle: { fontSize: 9 }, itemGap: 12 },
      grid: { left: 40, right: 10, top: 30, bottom: 35 },
      xAxis: { type: 'category', data: lag_times.map(String), name: lagAxisLabel, nameLocation: 'middle', nameGap: 25, nameTextStyle: { fontSize: 11, fontWeight: 500 }, axisLabel: { fontSize: 9, interval: Math.max(0, Math.floor(lag_times.length / 8) - 1) } },
      yAxis: { type: 'value', name: 'S(t)', nameLocation: 'middle', nameGap: 35, nameTextStyle: { fontSize: 11, fontWeight: 500 }, min: 0, max: 1.05, axisLabel: { fontSize: 9 } },
      series,
    };
  }, [data, interactionData]);

  if (!option) {
    return (
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', minHeight: { xs: 140, sm: 180 } }}>
        <Typography color="text.secondary">No contact events to analyze</Typography>
      </Box>
    );
  }

  return <ReactECharts option={option} style={{ height: '100%', minHeight: 140, width: '100%' }} opts={{ renderer: 'canvas' }} notMerge={true} />;
});

// Residence Distribution
interface ResidenceDistributionProps {
  data: KineticsData | null;
  interactionData: InteractionData | null;
}

export const ResidenceDistribution = memo(function ResidenceDistribution({
  data,
  interactionData,
}: ResidenceDistributionProps) {
  const option = useMemo(() => {
    if (!data || !data.residence_distribution) return null;

    const { bins, counts } = data.residence_distribution;
    const units = interactionData?.params?.units || 'ns';
    const normalizeBy = interactionData?.params?.normalize_by || 'actual_time';
    const durationAxisLabel = normalizeBy === 'counts' ? 'Duration (fr)' : `Duration (${getUnitLabel(units)})`;
    const durationUnit = normalizeBy === 'counts' ? 'fr' : getUnitLabel(units);

    if (bins.length === 0 || counts.every(c => c === 0)) return null;

    const maxCount = Math.max(...counts, 1);
    const maxBin = Math.max(...bins, 1);

    return {
      tooltip: { trigger: 'axis', formatter: (params: EChartsTooltipParam | EChartsTooltipParam[]) => {
        const paramsArray = Array.isArray(params) ? params : [params];
        if (paramsArray.length === 0) return '';
        const duration = paramsArray[0].axisValue;
        const count = paramsArray[0].data;
        return `Duration: ${duration}${durationUnit ? ' ' + durationUnit : ''}<br/>Count: ${count}`;
      }},
      grid: { left: 40, right: 10, top: 10, bottom: 35 },
      xAxis: { type: 'category', data: bins.map(String), name: durationAxisLabel, nameLocation: 'middle', nameGap: 25, nameTextStyle: { fontSize: 11, fontWeight: 500 }, axisLabel: { fontSize: 9, interval: Math.max(0, Math.ceil(maxBin / 10) - 1) } },
      yAxis: { type: 'value', name: 'Count', nameLocation: 'middle', nameGap: 35, nameTextStyle: { fontSize: 11, fontWeight: 500 }, axisLabel: { fontSize: 9 }, minInterval: 1, splitNumber: Math.min(maxCount, 5) },
      series: [{ type: 'bar', data: counts, itemStyle: { color: colors.data.query, borderColor: colors.data.queryDark, borderWidth: 1 }, barGap: '10%' }],
    };
  }, [data, interactionData]);

  if (!option) {
    return (
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', minHeight: { xs: 140, sm: 180 } }}>
        <Typography color="text.secondary">No residence events to display</Typography>
      </Box>
    );
  }

  return <ReactECharts option={option} style={{ height: '100%', minHeight: 140, width: '100%' }} opts={{ renderer: 'canvas' }} notMerge={true} />;
});
