/**
 * TypeScript types for ECharts components.
 * Provides type safety for chart callbacks and formatters.
 */

import type { ECharts } from 'echarts';

/** ECharts tooltip callback parameter for single data point */
export interface EChartsTooltipParam {
  componentType: string;
  seriesType: string;
  seriesIndex: number;
  seriesName: string;
  name: string;
  dataIndex: number;
  data: unknown;
  value: unknown;
  color: string;
  axisValue?: string | number;
  axisValueLabel?: string;
  marker?: string;
}

/** Heatmap data point: [x, y, value] */
export type HeatmapDataPoint = [number, number, number];

/** ECharts click event parameter */
export interface EChartsClickEventParams {
  componentType: string;
  seriesType: string;
  seriesIndex: number;
  seriesName: string;
  name: string;
  dataIndex: number;
  data: unknown;
  value: unknown;
  color?: string;
  event?: {
    offsetX: number;
    offsetY: number;
    event: MouseEvent;
  };
}

/** ZRender event parameter for low-level canvas interactions */
export interface ZRenderEventParams {
  offsetX: number;
  offsetY: number;
  event: MouseEvent;
  target?: unknown;
  topTarget?: unknown;
}

/** ECharts instance type */
export type EChartsInstance = ECharts;
