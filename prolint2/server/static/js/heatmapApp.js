///////////////////////////////////////////
///////////// Heatmap App /////////////////
///////////////////////////////////////////

export function heatmapApp(responseData) {

    var startLipidID = responseData['tableData'][0].lipidID;
    var startResidueID = responseData['topLipids'][0]

    var heatmapRoot = am5.Root.new("chartdiv4");

    // Set themes
    heatmapRoot.setThemes([
        am5themes_Animated.new(heatmapRoot)
    ]);

    // Create chart
    var heatmapChart = heatmapRoot.container.children.push(am5xy.XYChart.new(heatmapRoot, {
        panX: false,
        panY: false,
        wheelX: "none",
        wheelY: "none",
        layout: heatmapRoot.verticalLayout
    }));

    heatmapChart.children.unshift(am5.Label.new(heatmapRoot, {
        text: `Interactions between ResidueID: ${startResidueID} and LipidID: ${startLipidID}`,
        id: "besiTest",
        x: am5.p50,
        centerX: am5.percent(47),
        centerY: am5.percent(25),
        position: "absolute",
        fontWeight: "500",
        fontStyle: "oblique",
    }));

    // Create axes and their renderers
    var hmYRenderer = am5xy.AxisRendererY.new(heatmapRoot, {
        visible: false,
        minGridDistance: 20,
        inversed: true
    });

    hmYRenderer.grid.template.set("visible", false);

    var hmYAxis = heatmapChart.yAxes.push(am5xy.CategoryAxis.new(heatmapRoot, {
        maxDeviation: 0,
        renderer: hmYRenderer,
        categoryField: "ResidueAtoms"
    }));
    hmYAxis.children.moveValue(am5.Label.new(heatmapRoot, {
        text: "Residue atoms",
        rotation: -90,
        position: "absolute",
        y: am5.percent(80),
        // x: am5.percent(0),
        centerX: am5.percent(10),
        centerY: am5.percent(75),
        fontWeight: "500",
        // fontStyle: "oblique",
    }), 0);

    var hmXRenderer = am5xy.AxisRendererX.new(heatmapRoot, {
        visible: false,
        minGridDistance: 30,
        opposite: true
    });

    hmXRenderer.grid.template.set("visible", false);

    var hmXAxis = heatmapChart.xAxes.push(am5xy.CategoryAxis.new(heatmapRoot, {
        renderer: hmXRenderer,
        categoryField: "LipidAtoms"
    }));
    hmXAxis.children.moveValue(am5.Label.new(heatmapRoot, {
        text: "Lipid atoms",
        position: "absolute",
        // top x label :
        // x: am5.percent(50),
        // centerX: am5.percent(50),
        // centerY: am5.percent(70),

        // bottom x label
        x: am5.percent(50),
        y: am5.percent(1025),
        centerX: am5.percent(50),

        fontWeight: "500",
        // fontStyle: "oblique",
    }), hmXAxis.children.length - 1);

    // Create series
    var heatmapSeries = heatmapChart.series.push(am5xy.ColumnSeries.new(heatmapRoot, {
        calculateAggregates: true,
        stroke: am5.color(0xffffff),
        clustered: false,
        xAxis: hmXAxis,
        yAxis: hmYAxis,
        categoryXField: "LipidAtoms",
        categoryYField: "ResidueAtoms",
        valueField: "value"
    }));

    heatmapSeries.columns.template.setAll({
        tooltipText: "{value}",
        strokeOpacity: 1,
        strokeWidth: 2,
        width: am5.percent(100),
        height: am5.percent(100)
    });

    heatmapSeries.columns.template.events.on("pointerover", function (event) {
        var di = event.target.dataItem;
        if (di) {
            heatLegend.showValue(di.get("value", 0));
        }
    });

    heatmapSeries.events.on("datavalidated", function () {
        heatLegend.set("startValue", heatmapSeries.getPrivate("valueHigh"));
        heatLegend.set("endValue", heatmapSeries.getPrivate("valueLow"));
    });

    // Set up heat rules
    heatmapSeries.set("heatRules", [{
        target: heatmapSeries.columns.template,
        min: am5.color("#E74C3C"),
        max: am5.color("#FDEDEC"),
        dataField: "value",
        key: "fill"
    }]);

    // Add heat legend
    var heatLegend = heatmapChart.bottomAxesContainer.children.push(am5.HeatLegend.new(heatmapRoot, {
        orientation: "horizontal",
        startColor: am5.color("#FDEDEC"),
        endColor: am5.color("#E74C3C"),
    }));

    // Set data
    heatmapSeries.data.setAll(responseData['heatmapData']);
    hmYAxis.data.setAll(responseData['residueAtomsData']);
    hmXAxis.data.setAll(responseData['lipidAtomsData']);

    // Make stuff animate on load
    heatmapChart.appear(1000, 100);

    return [heatmapRoot, heatmapChart, heatmapSeries, hmYAxis, hmXAxis, heatLegend]

}