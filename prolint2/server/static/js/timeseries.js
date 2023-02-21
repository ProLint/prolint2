export function timeSeriesApp(contactData) {
    // Create root element
    // https://www.amcharts.com/docs/v5/getting-started/#Root_element
    var root = am5.Root.new("chartdiv_z");


    // Set themes
    // https://www.amcharts.com/docs/v5/concepts/themes/
    root.setThemes([
    am5themes_Animated.new(root)
    ]);


    // Create chart
    // https://www.amcharts.com/docs/v5/charts/xy-chart/
    var chart = root.container.children.push(am5xy.XYChart.new(root, {
    panX: true,
    panY: true,
    wheelX: "panX",
    wheelY: "zoomX",
    // scrollbarX: am5.Scrollbar.new(root, { orientation: "horizontal" }),
    // scrollbarY: am5.Scrollbar.new(root, { orientation: "vertical" }),
    pinchZoomX: true
    }));


    // Add cursor
    // https://www.amcharts.com/docs/v5/charts/xy-chart/cursor/
    var cursor = chart.set("cursor", am5xy.XYCursor.new(root, {}));
    cursor.lineY.set("visible", false);


    // Create axes
    // https://www.amcharts.com/docs/v5/charts/xy-chart/axes/
    var xRenderer = am5xy.AxisRendererX.new(root, {
    minGridDistance: 15
    });

    xRenderer.labels.template.setAll({
    rotation: -90,
    centerY: am5.p50,
    centerX: 0
    });

    xRenderer.grid.template.setAll({
    visible: true
    });

    var xAxis = chart.xAxes.push(am5xy.CategoryAxis.new(root, {
    maxDeviation: 0.3,
    categoryField: "category",
    renderer: xRenderer,
    tooltip: am5.Tooltip.new(root, {})
    }));

    var yAxis = chart.yAxes.push(am5xy.ValueAxis.new(root, {
    maxDeviation: 0.3,
    renderer: am5xy.AxisRendererY.new(root, {})
    }));


    // Create series
    // https://www.amcharts.com/docs/v5/charts/xy-chart/series/
    var series = chart.series.push(am5xy.ColumnSeries.new(root, {
    xAxis: xAxis,
    yAxis: yAxis,
    valueYField: "value",
    categoryXField: "category",
    adjustBulletPosition: false,
    tooltip: am5.Tooltip.new(root, {
        labelText: "{valueY}"
    })
    }));
    series.columns.template.setAll({
    width: 0.5
    });

    series.bullets.push(function() {
    return am5.Bullet.new(root, {
        locationY: 1,
        sprite: am5.Circle.new(root, {
        radius: 5,
        fill: series.get("fill")
        })
    })
    })

    var data = getTimeData(contactData);

    xAxis.data.setAll(data);
    series.data.setAll(data);


    // Make stuff animate on load
    // https://www.amcharts.com/docs/v5/concepts/animations/
    series.appear(1000);
    chart.appear(1000, 100);

    return [xAxis, series]

}

export function getTimeData(contactData) {
    var names = [];
    var data = [];
    for (var i = 0; i < contactData.length; i++) {
        names.push(contactData[i].residue);
        data.push({ category: contactData[i].residue, value: contactData[i].value });
    }
    return data
}