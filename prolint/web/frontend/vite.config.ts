import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import path from 'path'

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react()],
  resolve: {
    alias: {
      // Allow importing from the config theme directory
      '@prolint/config': path.resolve(__dirname, '../../config'),
    },
  },
  build: {
    rollupOptions: {
      output: {
        // Split heavy dependencies into separate chunks for better caching
        manualChunks(id) {
          if (id.includes('node_modules')) {
            // Molstar and its dependencies
            if (id.includes('molstar') || id.includes('immer') || id.includes('rxjs')) {
              return 'molstar';
            }
            // ECharts
            if (id.includes('echarts')) {
              return 'echarts';
            }
            // Cytoscape
            if (id.includes('cytoscape')) {
              return 'cytoscape';
            }
            // MUI and Emotion
            if (id.includes('@mui') || id.includes('@emotion')) {
              return 'mui';
            }
          }
        },
      },
    },
    chunkSizeWarningLimit: 1000,
  },
  server: {
    port: 3000,
    proxy: {
      '/api': {
        target: 'http://localhost:8000',
        changeOrigin: true,
      },
      '/ws': {
        target: 'ws://localhost:8000',
        ws: true,
      },
    },
    // Allow serving files from shared directory
    fs: {
      allow: ['..', '../..'],
    },
  },
})
