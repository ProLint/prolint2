/** ProLint2 Interactive Dashboard - Two Page Layout with lazy-loaded routes. */

import React, { Suspense } from 'react';
import { BrowserRouter as Router, Routes, Route, Navigate } from 'react-router-dom';
import { ThemeProvider, createTheme } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import Box from '@mui/material/Box';
import CircularProgress from '@mui/material/CircularProgress';

import DashboardLayout from './components/DashboardLayout';
import { colors } from './theme/visualizationTheme';

// Lazy load pages for code splitting - ViewPage has heavy deps (molstar, echarts, cytoscape)
const ComputePage = React.lazy(() => import('./pages/ComputePage'));
const ViewPage = React.lazy(() => import('./pages/ViewPage'));

/** Loading fallback for lazy-loaded routes. */
function PageLoader() {
  return (
    <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '50vh' }}>
      <CircularProgress />
    </Box>
  );
}

// ProLint unified theme - Modern scientific color scheme
const theme = createTheme({
  palette: {
    mode: 'light',
    primary: {
      light: colors.primary[400],
      main: colors.primary[500],
      dark: colors.primary[700],
    },
    secondary: {
      light: colors.accent[400],
      main: colors.accent[500],
      dark: colors.accent[600],
    },
    success: {
      light: colors.success.light,
      main: colors.success.main,
      dark: colors.success.dark,
    },
    warning: {
      light: colors.warning.light,
      main: colors.warning.main,
      dark: colors.warning.dark,
    },
    error: {
      light: colors.error.light,
      main: colors.error.main,
      dark: colors.error.dark,
    },
    info: {
      light: colors.info.light,
      main: colors.info.main,
      dark: colors.info.dark,
    },
    background: {
      default: colors.background.default,
      paper: colors.background.paper,
    },
    text: {
      primary: colors.text.primary,
      secondary: colors.text.secondary,
      disabled: colors.text.disabled,
    },
    divider: colors.border.default,
  },
  typography: {
    fontFamily: [
      '-apple-system',
      'BlinkMacSystemFont',
      '"Segoe UI"',
      'Roboto',
      '"Helvetica Neue"',
      'Arial',
      'sans-serif',
    ].join(','),
  },
  shape: {
    borderRadius: 8,
  },
});

function App() {
  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <Router>
        <Box sx={{ display: 'flex', minHeight: '100vh' }}>
          <DashboardLayout>
            <Suspense fallback={<PageLoader />}>
              <Routes>
                <Route path="/" element={<Navigate to="/compute" replace />} />
                <Route path="/compute" element={<ComputePage />} />
                <Route path="/view" element={<ViewPage />} />
              </Routes>
            </Suspense>
          </DashboardLayout>
        </Box>
      </Router>
    </ThemeProvider>
  );
}

export default App;
