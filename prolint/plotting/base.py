"""Plotting base module.

This module provides the base classes and registry for ProLint plotters.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Type
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from prolint.analysis.base import AnalysisResult


class BasePlotter(ABC):
    """Abstract base class for all ProLint plotters.

    Plotters convert AnalysisResult objects into matplotlib visualizations.
    Subclasses must implement the ``plot()`` method.

    Attributes
    ----------
    name : str
        Plotter name for registry lookup.
    required_analysis : str
        Name of the analysis type this plotter expects.
    description : str
        Human-readable description.

    See Also
    --------
    PlottingRegistry : Registry for creating plotters by name
    plot : Convenience function for plotting
    """

    name: str = "base_plotter"
    required_analysis: str = ""
    description: str = "Base plotter class"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that the AnalysisResult contains required data keys.

        Parameters
        ----------
        result : AnalysisResult
            Result to validate.

        Raises
        ------
        ValueError
            If result is missing required keys for this plotter.
        """
        pass

    @classmethod
    @abstractmethod
    def plot(cls, result: AnalysisResult, **kwargs) -> Tuple[Figure, Axes]:
        """Create the plot from an AnalysisResult.

        Parameters
        ----------
        result : AnalysisResult
            Analysis result containing data to plot.
        **kwargs : dict
            Plotter-specific options.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        pass


class PlottingRegistry:
    """Registry for plotter types.

    Manages registration and creation of plotter classes. All built-in
    plotters are registered automatically on import.

    Examples
    --------
    List available plotters:

    >>> from prolint.plotting import PlottingRegistry
    >>> PlottingRegistry.available()
    ['heatmap', 'density_map', 'survival_curve', ...]

    Create a plot:

    >>> fig, ax = PlottingRegistry.plot("heatmap", result)
    """

    _registry: Dict[str, Type[BasePlotter]] = {}

    @classmethod
    def register(cls, name: str, plotter_class: Type[BasePlotter]) -> None:
        """Register a plotter class.

        Parameters
        ----------
        name : str
            Name to register under.
        plotter_class : type
            Plotter class (must inherit from BasePlotter).
        """
        if not issubclass(plotter_class, BasePlotter):
            raise TypeError(f"{plotter_class} must be a subclass of BasePlotter")
        cls._registry[name] = plotter_class

    @classmethod
    def plot(cls, name: str, result: AnalysisResult, **kwargs) -> Tuple[Figure, Axes]:
        """Create a plot using a registered plotter.

        Parameters
        ----------
        name : str
            Plotter type name.
        result : AnalysisResult
            Analysis result to plot.
        **kwargs : dict
            Plotter-specific options.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.

        Raises
        ------
        ValueError
            If plotter name is not registered.
        """
        if name not in cls._registry:
            available = ", ".join(sorted(cls._registry.keys()))
            raise ValueError(f"Unknown plot type: '{name}'. Available: {available}")

        plotter = cls._registry[name]
        plotter.validate_result(result)
        return plotter.plot(result, **kwargs)

    @classmethod
    def available(cls) -> List[str]:
        """List all available plot types.

        Returns
        -------
        list of str
            Registered plotter names.
        """
        return sorted(cls._registry.keys())

    @classmethod
    def get_class(cls, name: str) -> Type[BasePlotter]:
        """Get a plotter class by name.

        Parameters
        ----------
        name : str
            Plotter name.

        Returns
        -------
        type
            Plotter class.
        """
        if name not in cls._registry:
            raise ValueError(f"Unknown plot type: {name}")
        return cls._registry[name]

    @classmethod
    def get_info(cls, name: str) -> Dict[str, str]:
        """Get information about a plotter.

        Parameters
        ----------
        name : str
            Plotter name.

        Returns
        -------
        dict
            Dict with name, required_analysis, and description keys.
        """
        plotter = cls.get_class(name)
        return {
            "name": plotter.name,
            "required_analysis": plotter.required_analysis,
            "description": plotter.description,
        }

    @classmethod
    def list_plotters(cls) -> List[Dict[str, str]]:
        """List all plotters with their info.

        Returns
        -------
        list of dict
            Info dict for each registered plotter.
        """
        return [cls.get_info(name) for name in cls.available()]


def plot(name: str, result: AnalysisResult, **kwargs) -> Tuple[Figure, Axes]:
    """Create a plot from an AnalysisResult.

    Convenience function that delegates to PlottingRegistry.plot().

    Parameters
    ----------
    name : str
        Plotter type name.
    result : AnalysisResult
        Analysis result to plot.
    **kwargs : dict
        Plotter-specific options.

    Returns
    -------
    tuple of (Figure, Axes)
        Matplotlib figure and axes objects.

    Examples
    --------
    >>> from prolint.plotting import plot
    >>> fig, ax = plot("heatmap", timeseries_result)
    >>> fig, ax = plot("survival_curve", kinetics_result)
    """
    return PlottingRegistry.plot(name, result, **kwargs)
