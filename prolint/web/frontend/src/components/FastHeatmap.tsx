/** Fast ECharts-based heatmap with Canvas rendering for large datasets. */

import { useMemo, useRef, useEffect } from 'react';
import ReactECharts from 'echarts-for-react';
import { Box, Typography } from '@mui/material';
import { colorScales } from '../theme/visualizationTheme';
import type { EChartsInstance, EChartsClickEventParams, HeatmapDataPoint } from '../types/echarts';

interface FastHeatmapProps {
  data: number[][];
  rowLabels: string[];
  colLabels: string[];
  colorScheme?: 'viridis' | 'blues' | 'prolint' | 'mako';
  maxHeight?: number;
  height?: string | number;
  showRowLabels?: boolean;
  showColorbar?: boolean;
  showXAxis?: boolean;
  xAxisLabel?: string;
  yAxisLabel?: string;
  onRowClick?: (rowLabel: string, rowIndex: number) => void;
  disableSampling?: boolean;
  rowTotals?: number[];
}

// Color schemes - using unified theme colors
const COLOR_SCHEMES = {
  viridis: colorScales.viridis,
  blues: colorScales.blues,
  prolint: colorScales.prolint,
  mako: colorScales.mako,
};

export default function FastHeatmap({
  data,
  rowLabels,
  colLabels,
  colorScheme = 'viridis',
  maxHeight = 500,
  height,
  showRowLabels = true,
  showColorbar = true,
  showXAxis = true,
  xAxisLabel = '',
  yAxisLabel = '',
  onRowClick,
  disableSampling = false,
  rowTotals,
}: FastHeatmapProps) {

  const { chartData, minValue, maxValue, sampledColLabels, sampleStep } = useMemo(() => {
    // Sample columns if too many (unless disabled)
    const maxCols = 200;
    let step = 1;
    let sampledLabels = colLabels;
    let sampledData = data;

    if (!disableSampling && colLabels.length > maxCols) {
      step = Math.ceil(colLabels.length / maxCols);
      sampledLabels = colLabels.filter((_, i) => i % step === 0);
      sampledData = data.map(row => row.filter((_, i) => i % step === 0));
    }

    // Convert to ECharts format: [colIndex, rowIndex, value]
    const echartsData: [number, number, number][] = [];
    let min = Infinity;
    let max = -Infinity;

    sampledData.forEach((row, rowIdx) => {
      row.forEach((value, colIdx) => {
        echartsData.push([colIdx, rowIdx, value]);
        if (value < min) min = value;
        if (value > max) max = value;
      });
    });

    return {
      chartData: echartsData,
      minValue: min === Infinity ? 0 : min,
      maxValue: max === -Infinity ? 1 : max,
      sampledColLabels: sampledLabels,
      sampleStep: step,
    };
  }, [data, colLabels, disableSampling]);

  const numRows = rowLabels.length;
  const numCols = sampledColLabels.length;

  // Calculate height - use explicit height if provided, otherwise calculate based on rows
  const rowHeight = numRows <= 30 ? 18 : numRows <= 60 ? 12 : 8;
  const calculatedHeight = Math.min(Math.max(numRows * rowHeight + 100, 300), maxHeight);
  const plotHeight = height !== undefined ? height : calculatedHeight;

  if (numRows === 0 || numCols === 0) {
    return <Typography color="text.secondary">No data to display</Typography>;
  }

  const option = {
    tooltip: {
      position: 'top',
      formatter: (params: EChartsClickEventParams) => {
        const [colIdx, rowIdx, value] = params.data as HeatmapDataPoint;
        let tooltip = `<b>Residue ${rowLabels[rowIdx]}</b><br/>Frame: ${sampledColLabels[colIdx]}`;
        if (rowTotals && rowTotals[rowIdx] !== undefined) {
          tooltip += `<br/>Total Contacts: ${rowTotals[rowIdx]}`;
        }
        return tooltip;
      },
    },
    grid: {
      top: 10,
      right: 10,
      bottom: showColorbar ? 60 : (showXAxis ? 30 : 10),
      left: 15,
      containLabel: false,
    },
    xAxis: {
      type: 'category',
      data: sampledColLabels,
      name: showXAxis ? xAxisLabel : '',
      nameLocation: 'middle',
      nameGap: 20,
      nameTextStyle: {
        fontSize: 11,
        fontWeight: 500,
      },
      splitArea: { show: false },
      axisLabel: {
        show: showXAxis,
        interval: Math.max(0, Math.floor(numCols / 10) - 1),
        fontSize: 11,
        rotate: 0,
      },
      axisTick: { show: showXAxis },
      axisLine: { show: true },
    },
    yAxis: {
      type: 'category',
      data: rowLabels,
      name: yAxisLabel,
      nameLocation: 'middle',
      nameGap: showRowLabels && numRows <= 40 ? 10 : 5,
      nameTextStyle: {
        fontSize: 11,
        fontWeight: 500,
      },
      splitArea: { show: false },
      axisLabel: {
        show: showRowLabels && numRows <= 40,
        fontSize: 11,
        interval: 0,
      },
      axisTick: { show: showRowLabels && numRows <= 40 },
      axisLine: { show: true },
      inverse: true,
    },
    visualMap: {
      show: showColorbar,
      min: minValue,
      max: maxValue,
      calculable: true,
      orient: 'horizontal',
      left: 'center',
      bottom: -8,
      itemWidth: 10,
      itemHeight: 150,
      textStyle: { fontSize: 9 },
      inRange: {
        color: COLOR_SCHEMES[colorScheme],
      },
    },
    series: [
      {
        name: 'Contacts',
        type: 'heatmap',
        data: chartData,
        emphasis: {
          itemStyle: {
            shadowBlur: 10,
            shadowColor: 'rgba(0, 0, 0, 0.5)',
          },
        },
        progressive: 0,
        animation: false,
      },
    ],
  };

  // Store current onRowClick in a ref to avoid stale closures
  const onRowClickRef = useRef(onRowClick);
  const rowLabelsRef = useRef(rowLabels);

  useEffect(() => {
    onRowClickRef.current = onRowClick;
    rowLabelsRef.current = rowLabels;
  }, [onRowClick, rowLabels]);

  // Callback when chart is ready - bind click event
  const onChartReady = (chart: EChartsInstance) => {
    if (!chart) return;

    // Remove existing handler
    chart.off('click');

    // Add click handler
    chart.on('click', (params) => {
      if (onRowClickRef.current) {
        const clickParams = params as EChartsClickEventParams;
        const data = clickParams.data || clickParams.value;
        if (data && Array.isArray(data)) {
          const rowIdx = data[1] as number;
          const rowLabel = rowLabelsRef.current[rowIdx];
          onRowClickRef.current(rowLabel, rowIdx);
        }
      }
    });
  };

  const containerStyle = height === '100%'
    ? { height: '100%', display: 'flex', flexDirection: 'column' as const }
    : {};

  // Generate accessible description
  const accessibleDescription = `Heatmap showing ${numRows} ${yAxisLabel || 'rows'} by ${colLabels.length} ${xAxisLabel || 'columns'}. Values range from ${minValue.toFixed(2)} to ${maxValue.toFixed(2)}.`;

  return (
    <Box
      sx={containerStyle}
      role="img"
      aria-label={accessibleDescription}
    >
      <ReactECharts
        option={option}
        style={{ height: plotHeight, width: '100%', flex: height === '100%' ? 1 : undefined }}
        opts={{ renderer: 'canvas' }}
        notMerge={true}
        onChartReady={onChartReady}
      />
      {sampleStep > 1 && (
        <Typography variant="caption" color="text.secondary" sx={{ display: 'block', textAlign: 'center', mt: 0.5 }}>
          Showing every {sampleStep}th frame ({numCols} of {colLabels.length} frames)
        </Typography>
      )}
    </Box>
  );
}
