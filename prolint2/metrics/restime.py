import numpy as np
try:
    from scipy.optimize import curve_fit
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available. Garcia & Stiller method will fall back to simple mean.")

from prolint2.metrics.base import FittingFunction
from prolint2.metrics.fitters import FittingFunctionFactory


class GarciaStillerSurvivalFunction:
    """
    Garcia & Stiller (1993) survival time correlation function implementation.
    This is the scientific standard for MD simulations residence time analysis.
    
    Formula: σ(t) = (1/N_j) * (1/(T-t)) * Σ[excess times] / σ(0)
    """
    
    def __init__(self, durations, t_total, delta_t_list):
        self.durations = durations
        self.t_total = t_total
        self.delta_t_list = delta_t_list
        self.survival_func = self._calculate_survival_function()
    
    def _calculate_survival_function(self):
        num_of_contacts = len(self.durations)
        if num_of_contacts == 0:
            return {delta_t: 0 for delta_t in self.delta_t_list}
            
        survival_func = {}
        
        for delta_t in self.delta_t_list:
            if delta_t == 0:
                # At t=0, survival function equals 1
                survival_func[delta_t] = 1.0
                # Calculate σ(0) normalization factor
                survival_func0 = float(sum([res_time - delta_t for res_time in self.durations if res_time >= delta_t])) / \
                         ((self.t_total - delta_t) * num_of_contacts)
            else:
                try:
                    # For t > 0, normalize by σ(0)
                    if survival_func0 > 0:
                        survival_func[delta_t] = float(sum([res_time - delta_t for res_time in self.durations if res_time >= delta_t])) / \
                                         ((self.t_total - delta_t) * num_of_contacts * survival_func0)
                    else:
                        survival_func[delta_t] = 0.0
                except (ZeroDivisionError, UnboundLocalError):
                    survival_func[delta_t] = 0.0
        
        return survival_func


class SurvivalFunction:
    """
    Legacy survival function implementation.
    Kept for compatibility - use GarciaStillerSurvivalFunction for new analyses.
    """
    def __init__(self, durations, t_total, delta_t_list):
        self.durations = durations
        self.t_total = t_total
        self.delta_t_list = delta_t_list
        self.num_of_contacts = len(durations)
        self.survival_func = self.calculate()

    def _calc_survival_value(self, delta_t):
        """
        Calculate survival function S(t) = P(T > t)
        Returns the probability that a contact duration is greater than delta_t
        """
        if self.num_of_contacts == 0:
            return 0.0
        
        # Count contacts that last LONGER than delta_t (not >=)
        surviving_count = sum(1 for duration in self.durations if duration > delta_t)
        
        return surviving_count / self.num_of_contacts

    def calculate(self):
        survival_func = {}
        for delta_t in self.delta_t_list:
            survival_value = self._calc_survival_value(delta_t)
            survival_func[delta_t] = survival_value

        return survival_func


class KoffCalculator:
    """
    Garcia & Stiller (1993) based koff calculator - now the default method.
    Uses Garcia & Stiller survival time correlation function with bi-exponential fitting.
    """
    def __init__(
        self, durations, t_total, timestep, fitting_func_name="survival_correlation", 
        initial_guess=[1., 1., 1., 1.], cap=True, **kwargs
    ):
        self.durations = durations
        self.t_total = t_total
        self.timestep = timestep
        self.initial_guess = initial_guess
        self.cap = cap
        self.kwargs = kwargs
        self.fitting_func_name = fitting_func_name

        if self._is_empty_or_zeros(self.durations):
            self.res_time, self.koff = 0, 0
            return

        # Use Garcia & Stiller method by default for scientific accuracy
        if fitting_func_name == "survival_correlation":
            self._calculate_garcia_stiller_method()
        else:
            # Legacy method for compatibility
            self._calculate_legacy_method(fitting_func_name)

    def _is_empty_or_zeros(self, array):
        return len(array) == 0 or np.all(array == 0)

    def _calculate_garcia_stiller_method(self):
        """
        Calculate koff using Garcia & Stiller (1993) method.
        This is the recommended approach (default).
        """
        # Create time points
        self.delta_t_list = np.arange(0, self.t_total, self.timestep)
        
        # Calculate survival function using Garcia & Stiller method
        survival_calc = GarciaStillerSurvivalFunction(self.durations, self.t_total, self.delta_t_list)
        self.survival_func = survival_calc.survival_func
        
        # Perform bi-exponential fitting
        survival_rates = np.array([self.survival_func[delta_t] for delta_t in self.delta_t_list])
        
        try:
            res_time, koff, r_squared = self._garcia_stiller_curve_fitting(survival_rates)
            
            # Apply cap if requested
            if self.cap and res_time > self.t_total:
                res_time = self.t_total
                koff = 1.0 / res_time if res_time > 0 else 0
            
            self.res_time = res_time
            self.koff = koff
            self.r_squared = r_squared
            
        except Exception as e:
            # Fallback to simple mean if fitting fails
            print(f"Garcia & Stiller fitting failed: {e}, using simple mean")
            mean_duration = np.mean(self.durations)
            self.res_time = mean_duration
            self.koff = 1.0 / mean_duration
            self.r_squared = 0.0

    @staticmethod
    def _bi_expo(x, k1, k2, A, B):
        """Bi-exponential function for fitting"""
        try:
            exp1 = np.exp(np.clip(-k1 * x, -700, 700))  # Prevent overflow
            exp2 = np.exp(np.clip(-k2 * x, -700, 700))
            return A * exp1 + B * exp2
        except:
            return np.zeros_like(x)

    def _garcia_stiller_curve_fitting(self, survival_rates):
        """Garcia & Stiller bi-exponential curve fitting"""
        if not SCIPY_AVAILABLE:
            raise RuntimeError("scipy is required for Garcia & Stiller method")
        
        # Clean data
        survival_rates = np.nan_to_num(survival_rates)
        
        try:
            popt, pcov = curve_fit(
                self._bi_expo, 
                np.array(self.delta_t_list), 
                np.array(survival_rates), 
                p0=self.initial_guess, 
                maxfev=100000
            )
            
            # Calculate R-squared
            n_fitted = self._bi_expo(np.array(self.delta_t_list), *popt)
            ss_res = np.sum((np.nan_to_num(n_fitted) - np.nan_to_num(survival_rates))**2)
            ss_tot = np.sum((survival_rates - np.mean(survival_rates))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            # Extract koff (use smaller rate constant)
            ks = [abs(k) for k in popt[:2]]
            ks.sort()
            koff = ks[0] if ks[0] > 0 else 0
            res_time = 1/koff if koff > 0 else 0
            
            return res_time, koff, r_squared
            
        except Exception as e:
            raise RuntimeError(f"Curve fitting failed: {e}")

    def _calculate_legacy_method(self, fitting_func_name):
        """
        Legacy method using the original ProLint2 approach.
        Kept for compatibility.
        """
        if fitting_func_name not in FittingFunction.registry:
            func_names = ", ".join(FittingFunction.registry.keys())
            raise ValueError(
                f"Invalid fitting_func_name: {fitting_func_name}. Valid names are: {func_names}"
            )

        self.fitting_func = FittingFunctionFactory.get_fitting_function(fitting_func_name)
        
        # Improve time range for fitting
        self.delta_t_list = self._create_optimal_time_range()
        
        # Use legacy survival function
        survival_calc = SurvivalFunction(
            self.durations, np.max(self.t_total), self.delta_t_list
        )
        self.survival_func = survival_calc.survival_func
        self.res_time, self.koff = self._calculate_koff_legacy()

    def _create_optimal_time_range(self):
        """
        Create an optimal time range for fitting that avoids numerical issues.
        Use only the meaningful range where survival function changes.
        """
        if len(self.durations) == 0:
            return np.array([0])
        
        max_duration = max(self.durations)
        
        # Use time range from 0 to 2-3 times the maximum duration
        # This ensures we capture the decay without too many zeros
        max_fit_time = min(self.t_total, max_duration * 3)
        
        # Use reasonable number of points for fitting (50-200)
        # Avoid thousands of points that can cause numerical issues
        optimal_timestep = max(self.timestep, max_fit_time / 100)
        
        return np.arange(0, max_fit_time, optimal_timestep)

    def _calculate_koff_legacy(self):
        """Legacy koff calculation method"""
        survival_rates = np.array([
            self.survival_func[delta_t] for delta_t in self.delta_t_list
        ])
        
        # Check for degenerate cases that would cause fitting issues
        if len(survival_rates) < 3:
            # Too few points for meaningful fitting
            return np.mean(self.durations), 1.0 / np.mean(self.durations)
        
        non_zero_count = np.sum(survival_rates > 0)
        if non_zero_count < len(survival_rates) * 0.3:
            # Too many zeros - use simple average instead
            mean_duration = np.mean(self.durations)
            return mean_duration, 1.0 / mean_duration
        
        # Remove trailing zeros for better fitting
        last_nonzero = np.where(survival_rates > 0)[0]
        if len(last_nonzero) > 0:
            cutoff = last_nonzero[-1] + 1
            survival_rates = survival_rates[:cutoff]
            delta_t_fit = self.delta_t_list[:cutoff]
        else:
            delta_t_fit = self.delta_t_list
        
        try:
            popt = self.fitting_func.fit(
                np.array(delta_t_fit), survival_rates, **self.kwargs
            )
            koff = self.fitting_func.get_koff(popt)
            
            # Sanity check: koff should be positive and reasonable
            if koff <= 0 or koff > 1000:  # Arbitrary upper bound
                raise ValueError("Unreasonable koff value")
                
            res_time = 1 / koff
            
            # Another sanity check: residence time shouldn't be much larger than max duration
            max_duration = max(self.durations) if self.durations else 1
            if res_time > max_duration * 10:
                # Fall back to simple average
                mean_duration = np.mean(self.durations)
                return mean_duration, 1.0 / mean_duration
                
            return res_time, koff
            
        except Exception:
            # If fitting fails, fall back to simple average
            mean_duration = np.mean(self.durations)
            return mean_duration, 1.0 / mean_duration

    def calculate_koff(self):
        """
        Public method for backward compatibility.
        Returns the calculated residence time and koff.
        """
        return self.res_time, self.koff


# Backward compatibility - alias for users who might import the old class
class LegacyKoffCalculator(KoffCalculator):
    """
    Backward compatibility class for legacy ProLint2 method.
    Use: LegacyKoffCalculator(..., fitting_func_name="bi_expo")
    """
    def __init__(self, durations, t_total, timestep, fitting_func_name="bi_expo", **kwargs):
        super().__init__(durations, t_total, timestep, fitting_func_name, **kwargs)


# Helper function for easy method selection
def calculate_residence_time(durations, t_total, timestep, method="survival_correlation", **kwargs):
    """
    Calculate residence time using different methods.
    
    Parameters:
    -----------
    durations : list or array
        Contact durations
    t_total : float
        Total simulation time
    timestep : float
        Time step for analysis
    method : str
        Method to use: "survival_correlation" (default), "bi_expo", "mono_expo"
    **kwargs : additional arguments
    
    Returns:
    --------
    tuple : (residence_time, koff)
    """
    calculator = KoffCalculator(durations, t_total, timestep, fitting_func_name=method, **kwargs)
    return calculator.res_time, calculator.koff
