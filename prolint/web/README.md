# ProLint Web Dashboard

Interactive web dashboard for visualizing and analyzing biomolecular contacts from molecular dynamics simulations.

## Tech Stack

**Backend:**
- FastAPI
- Uvicorn
- WebSockets

**Frontend:**
- React 18
- TypeScript
- Vite
- MUI (Material-UI)

**Visualization:**
- ECharts
- Recharts
- Cytoscape (network graphs)
- Mol* (3D molecular structures)

## Project Structure

```
web/
├── backend/
│   ├── api/           # REST API endpoints
│   ├── models/        # Data models
│   ├── services/      # Business logic
│   ├── websocket/     # WebSocket handlers
│   └── main.py        # FastAPI application
├── frontend/
│   ├── src/
│   │   ├── components/  # React components
│   │   ├── pages/       # Page components
│   │   ├── hooks/       # Custom React hooks
│   │   ├── types/       # TypeScript types
│   │   ├── utils/       # Utility functions
│   │   └── App.tsx      # Main application
│   ├── package.json
│   └── vite.config.ts
└── requirements.txt     # Python dependencies
```

## Prerequisites

- Python 3.8+
- Node.js 18+
- npm

## Installation

**Backend:**

```bash
cd prolint/web
pip install -r requirements.txt
```

**Frontend:**

```bash
cd prolint/web/frontend
npm install
```

## Running Locally

**Start Backend (port 8000):**

```bash
cd prolint/web/backend
python -m uvicorn prolint.web.backend.main:app --reload --host 0.0.0.0 --port 8000
```

Or run directly:

```bash
python prolint/web/backend/main.py
```

**Start Frontend (port 3000):**

```bash
cd prolint/web/frontend
npm run dev
```

Access the dashboard at `http://localhost:3000`

## API Documentation

Once the backend is running:

- **Swagger UI:** http://localhost:8000/api/docs
- **ReDoc:** http://localhost:8000/api/redoc
- **Health Check:** http://localhost:8000/health

## Key Features

- **Dataset Upload:** Upload topology and trajectory files for analysis
- **Contact Computation:** Compute molecular contacts with real-time progress updates
- **Interactive Heatmaps:** Visualize contact frequencies and patterns
- **Timeseries Analysis:** Track contact dynamics over simulation time
- **Network Graphs:** Explore shared contacts between molecules
- **3D Visualization:** View molecular structures with Mol*