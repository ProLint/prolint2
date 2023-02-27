import { networkApp } from "../js/network.js";
///////////////////////////////////////////
/////////////// RadarApp //////////////////
///////////////////////////////////////////
export function radarApp() {
    
    var root = am5.Root.new("chartdiv_x");

    const theme = am5.Theme.new(root);
    theme.rule("Label").set("fontSize", 10);
    theme.rule("Grid").set("strokeOpacity", 0.06);

    root.setThemes([
        am5themes_Animated.new(root),
        theme
    ]);

    root.fps = 30;

    // Params
    var innerRadius = 50;

    // Create chart
    var chart = root.container.children.push(am5radar.RadarChart.new(root, {
        panX: false,
        panY: false,
        wheelX: "none",
        wheelY: "none",
        innerRadius: am5.percent(innerRadius),
        radius: am5.percent(85),
        startAngle: 270 - 170,
        endAngle: 270 + 170
        // Right handed half circle:
        // startAngle: -90,
        // endAngle: 90
    }));

    // Add cursor
    var cursor = chart.set("cursor", am5radar.RadarCursor.new(root, {
        behavior: "zoomX",
        radius: am5.percent(innerRadius),
        innerRadius: -15
    }));
    cursor.lineY.set("visible", true);
    cursor.lineY.set("opacity", 0.5);

    // Create axes and their renderers
    var xRenderer = am5radar.AxisRendererCircular.new(root, {
        minGridDistance: 10
    });

    xRenderer.labels.template.setAll({
        radius: 5,
        textType: "radial",
        centerY: am5.p50
    });

    // Metric axis
    var yRenderer = am5radar.AxisRendererRadial.new(root, {
        axisAngle: 90
    });
    // Metric axis labels
    yRenderer.labels.template.setAll({
        centerX: am5.p50
    });

    var categoryAxis = chart.xAxes.push(am5xy.CategoryAxis.new(root, {
        maxDeviation: 0,
        categoryField: "residue",
        renderer: xRenderer
    }));

    var modGridCategoryAxis = categoryAxis.get("renderer");
    modGridCategoryAxis.grid.template.setAll({
        strokeDasharray: [2, 2]
    });

    var valueAxis = chart.yAxes.push(am5xy.ValueAxis.new(root, {
        min: 0,
        strictMinMax: true,
        extraMax: 0.1,
        renderer: yRenderer
    }));

    var modGridValueAxis = valueAxis.get("renderer");
    modGridValueAxis.grid.template.setAll({
        strokeDasharray: [2, 2]
    });

    // We need to get mouse selection events to form a selection
    // which we'll pass to the 3D viewer.
    var viewerResidueSelection = {
        "start": undefined,
        "end": undefined,
    }
    cursor.events.on("selectstarted", function (ev) {
        var x = ev.target.getPrivate("positionX");
        var residue_id = categoryAxis.axisPositionToIndex(categoryAxis.toAxisPosition(x));
        viewerResidueSelection["start"] = residue_id
    });

    cursor.events.on("selectended", function (ev) {
        var x = ev.target.getPrivate("positionX");
        var residue_id = categoryAxis.axisPositionToIndex(categoryAxis.toAxisPosition(x));
        viewerResidueSelection["end"] = residue_id

        const residueRange = [
            parseInt(viewerResidueSelection["start"]),
            parseInt(viewerResidueSelection["end"])
        ];
        residueRange.sort((a, b) => a - b);

        // Only ok for very small selections.
        // Very slow otherwise.
        // selSec = []
        // for (let ix = residueRange[0]; ix < residueRange[1]; ix++) {
        //     columnColor = series.columns.values[ix].get("fill")
        //     selSec.push({
        //         residue_number: ix,
        //         representation: 'spacefill',
        //         representationColor: columnColor,
        //     })
        // }

        var selectSections = [{
            start_residue_number: residueRange[0],
            end_residue_number: residueRange[1],
            color: {
                r: 255,
                g: 0,
                b: 255
            },
            representation: 'spacefill',
            focus: true,
        }]
        viewerInstance.visual.select({
            data: selectSections, // selSec is very slow
        })

    });

    chart.zoomOutButton.events.on('click', function (ev) {
        viewerInstance.visual.clearSelection()
        viewerInstance.visual.reset({
            camera: true
        })
    })

    // Pie chart label
    var axisRange = categoryAxis.createAxisRange(categoryAxis.makeDataItem({
        above: true
    }));

    // Create series
    var series = chart.series.push(am5radar.RadarColumnSeries.new(root, {
        calculateAggregates: true,
        name: "Series",
        xAxis: categoryAxis,
        yAxis: valueAxis,
        valueYField: "value",
        categoryXField: "residue",
        tooltip: am5.Tooltip.new(root, {
            labelText: "{categoryX}: {valueY}"
        })
    }));

    series.on("tooltipDataItem", function (tooltipDataItem) {
        radarSeriesLegend.showValue(tooltipDataItem.get("valueY"))

        var residueID = tooltipDataItem.dataContext.residue.split(' ')[1]
        residueID = parseInt(residueID)

        var selectSections = [{
            residue_number: residueID,
            color: {
                r: 255,
                g: 0,
                b: 255
            },
        }]
        viewerInstance.visual.highlight({
            data: selectSections,
        })
    })

    series.columns.template.set("strokeOpacity", 0);

    series.events.on("datavalidated", function () {
        radarSeriesLegend.set("endValue", series.getPrivate("valueYHigh"));
        radarSeriesLegend.set("startValue", series.getPrivate("valueYLow"));
    });

    // Set up heat rules
    series.set("heatRules", [{
        target: series.columns.template,
        key: "fill",
        min: am5.color(0x673AB7),
        max: am5.color(0xF44336),
        dataField: "valueY",
        key: "fill"
    }]);

    var radarSeriesLegend = chart.bottomAxesContainer.children.push(am5.HeatLegend.new(root, {
        orientation: "horizontal",
        startColor: am5.color(0x673AB7),
        endColor: am5.color(0xF44336),
        width: 160,
        x: am5.percent(37.5),
        // puts the legen inside the plot
        // position: "absolute",
        // x: am5.percent(37.5),
        // centerY: am5.percent(1350),
        opacity: 1,
        // startText: "",
        // endText: "",
    }));
    radarSeriesLegend.startLabel.setAll({
        fontSize: 12,
        fill: radarSeriesLegend.get("startColor"),
        isMeasured: false,
        paddingLeft: -10,
        paddingTop: -1,
    });
    radarSeriesLegend.endLabel.setAll({
        fontSize: 12,
        fill: radarSeriesLegend.get("endColor"),
        isMeasured: false,
        x: am5.percent(75),
        paddingTop: -1,
    });


    function createRange(name, lipidData, index) {
        axisRange.get("label").setAll({
            text: name
        });
        // first residue
        axisRange.set("category", lipidData[0].residue);
        // last residue
        axisRange.set("endCategory", lipidData[lipidData.length - 1].residue);

        // every 3rd color for a bigger contrast
        var fill = axisRange.get("axisFill");
        fill.setAll({
            toggleKey: "active",
            cursorOverStyle: "pointer",
            fill: am5.color(0x095256),
            // fill: colorSet.getIndex(3),
            visible: true,
            innerRadius: -15,
            cornerRadius: 15,
        });
        axisRange.get("grid").set("visible", false);

        var label = axisRange.get("label");
        label.setAll({
            fill: am5.color(0xffffff),
            textType: "circular",
            visible: true,
            radius: -12
        });

        fill.events.on("click", function (event) {
            var dataItem = event.target.dataItem;
            if (event.target.get("active")) {
                categoryAxis.zoom(0, 1);
            } else {
                categoryAxis.zoomToCategories(dataItem.get("category"), dataItem.get("endCategory"));
            }
        });
    }

    return {
        root: root,
        chart: chart,
        series: series,
        axisRange: axisRange,
        categoryAxis: categoryAxis,
        createRange: createRange,
        active: true
    }
}
