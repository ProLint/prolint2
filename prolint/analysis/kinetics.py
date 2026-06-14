"""Kinetics analysis for binding/unbinding dynamics."""

import logging
from typing import Optional, List, Dict, Literal
import numpy as np
import warnings

from prolint.analysis.base import BaseAnalysis, AnalysisResult

logger = logging.getLogger(__name__)


class KineticsAnalysis(BaseAnalysis):
    """Kinetics analysis for binding/unbinding dynamics.

    Computes binding kinetics metrics including on/off rates, residence times,
    and survival curves with optional exponential fits.

    Attributes
    ----------
    MIN_EVENTS_MONO : int
        Minimum events required for monoexponential fit (default: 5).
    MIN_EVENTS_BI : int
        Minimum events required for biexponential fit (default: 25).

    See Also
    --------
    TimeSeriesAnalysis : Contact counts over time
    """

    name = "kinetics"
    """Analysis name for registry."""

    description = "Binding kinetics, residence times, and survival curves"
    """Human-readable description."""

    MIN_EVENTS_MONO = 5
    MIN_EVENTS_BI = 25

    def run(
        self,
        query_residue: int,
        database_residue: Optional[int] = None,
        database_type: Optional[str] = None,
        mode: Literal["individual", "accumulated"] = "individual",
        fit_survival: bool = True,
        max_lag: int = 100,
    ) -> AnalysisResult:
        """Compute kinetics analysis for a query residue.

        Parameters
        ----------
        query_residue : int
            Query residue ID to analyze.
        database_residue : int, optional
            Specific database residue ID. Required for "individual" mode.
        database_type : str, optional
            Database residue name (e.g., "CHOL"). Required for "accumulated" mode.
        mode : {"individual", "accumulated"}, default="individual"
            Analysis mode:

            - "individual": Single residue-residue pair kinetics
            - "accumulated": Aggregated kinetics across all molecules of a type
        fit_survival : bool, default=True
            Whether to fit exponential models to survival curves.
        max_lag : int, default=100
            Maximum lag time for survival curve computation.

        Returns
        -------
        AnalysisResult
            Result with data containing:

            - mode : str analysis mode
            - kinetics : dict with koff, kon, kd, residence_times,
              occupancy, n_events, n_frames
            - survival_curve : dict with lag_times, survival_probability,
              mono_fit, bi_fit, selected_model
            - residence_distribution : dict with bins and counts
            - contact_frames : list of frame indices with contacts

        Raises
        ------
        ValueError
            If database_residue not provided for "individual" mode, or
            database_type not provided for "accumulated" mode.
        """
        if mode == "individual" and database_residue is None:
            raise ValueError("database_residue required for 'individual' mode")
        if mode == "accumulated" and database_type is None:
            raise ValueError("database_type required for 'accumulated' mode")

        logger.info(
            "Computing kinetics for residue %d (mode=%s)",
            query_residue,
            mode,
        )

        n_frames = self.universe.trajectory.n_frames
        norm_factor = self.contacts.norm_factor

        # Durations are already scaled by norm_factor (user's chosen units)
        durations = []
        contacts_data = self.contacts.contacts
        if query_residue in contacts_data:
            query_data = contacts_data[query_residue]
            if mode == "individual" and database_residue is not None:
                # Find the database_residue in the nested structure
                db_id_to_resname = self._get_database_id_to_resname()
                target_resname = db_id_to_resname.get(database_residue)
                if target_resname and target_resname in query_data:
                    if database_residue in query_data[target_resname]:
                        durations = [
                            float(d)
                            for d in query_data[target_resname][database_residue]
                            if d > 0
                        ]
            elif mode == "accumulated" and database_type is not None:
                if database_type in query_data:
                    for dur_array in query_data[database_type].values():
                        durations.extend(
                            [float(d) for d in dur_array if d > 0]
                        )

        # Calculate occupancy based on mode
        if mode == "individual" and database_residue is not None:
            # For individual mode: occupancy = frames with contact / total frames
            contact_frames_data = self.contacts.contact_frames.get(query_residue, {})
            contact_frames_set = set(contact_frames_data.get(database_residue, []))
            occupancy = len(contact_frames_set) / n_frames if n_frames > 0 else 0.0
        elif mode == "accumulated" and database_type is not None:
            # For accumulated mode: use compute_metric to get occupancy across all residues of that type
            occ_result = self.contacts.compute_metric(
                "occupancy", target_resname=database_type
            )
            occupancy = 0.0
            if (
                query_residue in occ_result
                and database_type in occ_result[query_residue]
            ):
                occupancy = occ_result[query_residue][database_type]["global"]
        else:
            occupancy = 0.0

        # Scale max_lag to the same units as durations
        scaled_max_lag = max_lag * norm_factor

        kinetics = self._compute_kinetics(durations, occupancy, n_frames)
        survival_curve = self._compute_survival_curve(
            durations, scaled_max_lag, fit_survival, kinetics["koff"]
        )
        residence_dist = self._compute_residence_distribution(durations)

        # Collect contact frames using base class helper or direct access
        all_frames: set = set()
        if mode == "individual" and database_residue is not None:
            frames_data = self.contacts.contact_frames.get(query_residue, {})
            all_frames = set(frames_data.get(database_residue, []))
        elif mode == "accumulated" and database_type is not None:
            filtered = self._filter_by_database_type(database_type)
            for frames in filtered.get(query_residue, {}).values():
                all_frames.update(frames)

        return AnalysisResult(
            data={
                "mode": mode,
                "kinetics": kinetics,
                "survival_curve": survival_curve,
                "residence_distribution": residence_dist,
                "contact_frames": sorted(all_frames),
            },
            metadata={
                "query_residue": query_residue,
                "database_residue": database_residue,
                "database_type": database_type,
                "n_frames": n_frames,
            },
        )

    def _compute_kinetics(
        self, durations: List[float], occupancy: float, n_frames: int
    ) -> Dict:
        """Compute kinetics metrics from event durations (in user's chosen units)."""
        n_events = len(durations)
        mean_residence = float(np.mean(durations)) if durations else 0.0
        std_residence = float(np.std(durations)) if durations else 0.0
        max_residence = float(max(durations)) if durations else 0.0

        koff = 1.0 / mean_residence if mean_residence > 0 else 0.0
        n_non_contact = int((1 - occupancy) * n_frames)
        kon = n_events / n_non_contact if n_non_contact > 0 else 0.0
        kd = koff / kon if kon > 0 else None

        return {
            "koff": float(koff),
            "kon": float(kon),
            "kd": float(kd) if kd is not None else None,
            "mean_residence_time": float(mean_residence),
            "std_residence_time": float(std_residence),
            "max_residence_time": float(max_residence),
            "occupancy": float(occupancy),
            "n_events": n_events,
            "n_frames": n_frames,
        }

    def _compute_survival_curve(
        self, durations: List[float], max_lag: float, fit: bool, koff_estimate: float
    ) -> Dict:
        """Compute survival curve and optional exponential fits."""
        n_events = len(durations)
        max_duration = max(durations) if durations else 1.0
        n_points = min(int(max_duration) + 1, int(max_lag))
        lag_times = np.linspace(0, min(max_duration, max_lag), max(n_points, 2)).tolist()

        survival_probability = [
            sum(1 for d in durations if d >= lag) / n_events if n_events > 0 else 0
            for lag in lag_times
        ]

        result = {
            "lag_times": lag_times,
            "survival_probability": survival_probability,
            "mono_fit": None,
            "bi_fit": None,
            "selected_model": None,
            "min_events_mono": self.MIN_EVENTS_MONO,
            "min_events_bi": self.MIN_EVENTS_BI,
        }

        if not fit:
            return result

        t_data = np.array(lag_times, dtype=float)
        s_data = np.array(survival_probability, dtype=float)
        fit_mask = t_data > 0
        t_fit, s_fit = t_data[fit_mask], s_data[fit_mask]

        if len(t_fit) >= 2 and n_events >= self.MIN_EVENTS_MONO:
            result["mono_fit"] = self._fit_monoexponential(
                t_fit, s_fit, t_data, koff_estimate
            )

        if len(t_fit) >= 4 and n_events >= self.MIN_EVENTS_BI:
            result["bi_fit"] = self._fit_biexponential(
                t_fit, s_fit, t_data, koff_estimate
            )

        result["selected_model"] = self._select_model(
            result["mono_fit"], result["bi_fit"]
        )
        return result

    def _fit_monoexponential(
        self,
        t_fit: np.ndarray,
        s_fit: np.ndarray,
        t_full: np.ndarray,
        koff_estimate: float,
    ) -> Optional[Dict]:
        """Fit monoexponential decay: S(t) = exp(-k * t)."""
        try:
            from scipy.optimize import curve_fit

            def monoexp(t, k):
                return np.exp(-k * t)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p0 = [koff_estimate if koff_estimate > 0 else 0.1]
                popt, _ = curve_fit(
                    monoexp, t_fit, s_fit, p0=p0, bounds=([0.0001], [10.0]), maxfev=1000
                )
                k_fitted = popt[0]

                s_pred = monoexp(t_fit, k_fitted)
                ss_res = np.sum((s_fit - s_pred) ** 2)
                ss_tot = np.sum((s_fit - np.mean(s_fit)) ** 2)
                r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

                n_points = len(t_fit)
                aic = (
                    2 * 1 + n_points * np.log(ss_res / n_points)
                    if ss_res > 0
                    else float("inf")
                )

                logger.debug(
                    "Monoexponential fit: k_off=%.4f, R²=%.3f",
                    k_fitted,
                    r2,
                )

                return {
                    "k_off": float(k_fitted),
                    "r_squared": float(r2),
                    "aic": float(aic),
                    "fitted_curve": monoexp(t_full, k_fitted).tolist(),
                    "half_life": float(np.log(2) / k_fitted) if k_fitted > 0 else None,
                }
        except Exception as e:
            logger.debug("Monoexponential fit failed: %s", e)
            return None

    def _fit_biexponential(
        self,
        t_fit: np.ndarray,
        s_fit: np.ndarray,
        t_full: np.ndarray,
        koff_estimate: float,
    ) -> Optional[Dict]:
        """Fit biexponential decay: S(t) = a*exp(-k1*t) + (1-a)*exp(-k2*t)."""
        try:
            from scipy.optimize import curve_fit

            def biexp(t, a, k1, k2):
                return a * np.exp(-k1 * t) + (1 - a) * np.exp(-k2 * t)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                k_fast = koff_estimate * 2 if koff_estimate > 0 else 0.2
                k_slow = koff_estimate * 0.5 if koff_estimate > 0 else 0.05
                p0 = [0.5, k_fast, k_slow]

                popt, _ = curve_fit(
                    biexp,
                    t_fit,
                    s_fit,
                    p0=p0,
                    bounds=([0.01, 0.001, 0.0001], [0.99, 10.0, 10.0]),
                    maxfev=2000,
                )
                a, k1, k2 = popt

                if k1 < k2:
                    a, k1, k2 = 1 - a, k2, k1

                s_pred = biexp(t_fit, a, k1, k2)
                ss_res = np.sum((s_fit - s_pred) ** 2)
                ss_tot = np.sum((s_fit - np.mean(s_fit)) ** 2)
                r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

                n_points = len(t_fit)
                aic = (
                    2 * 3 + n_points * np.log(ss_res / n_points)
                    if ss_res > 0
                    else float("inf")
                )

                logger.debug(
                    "Biexponential fit: k_fast=%.4f, k_slow=%.4f, R²=%.3f",
                    k1,
                    k2,
                    r2,
                )

                return {
                    "a_fast": float(a),
                    "k_fast": float(k1),
                    "k_slow": float(k2),
                    "r_squared": float(r2),
                    "aic": float(aic),
                    "fitted_curve": biexp(t_full, a, k1, k2).tolist(),
                    "half_life_fast": float(np.log(2) / k1) if k1 > 0 else None,
                    "half_life_slow": float(np.log(2) / k2) if k2 > 0 else None,
                }
        except Exception as e:
            logger.debug("Biexponential fit failed: %s", e)
            return None

    def _select_model(
        self, mono_fit: Optional[Dict], bi_fit: Optional[Dict]
    ) -> Optional[str]:
        """Select best model based on AIC."""
        if mono_fit and bi_fit:
            return (
                "biexponential"
                if mono_fit["aic"] - bi_fit["aic"] > 2
                else "monoexponential"
            )
        return "monoexponential" if mono_fit else ("biexponential" if bi_fit else None)

    def _compute_residence_distribution(
        self, durations: List[float], max_bins: int = 50
    ) -> Dict:
        """Compute histogram of residence time durations."""
        if not durations:
            return {"bins": [], "counts": []}

        n_bins = min(len(set(durations)), max_bins)
        counts_arr, bin_edges = np.histogram(durations, bins=max(n_bins, 1))
        # Use bin centers as labels
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        return {
            "bins": bin_centers.tolist(),
            "counts": counts_arr.tolist(),
        }
