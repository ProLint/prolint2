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
    panY: false,
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

    yAxis.children.moveValue(am5.Label.new(root, {
        rotation: -90,
        text: "Contact Metric",
        y: am5.p50,
        centerX: am5.p50
      }), 0);
      
    var modGridxAxis = xAxis.get("renderer");
    modGridxAxis.grid.template.setAll({
        strokeDasharray: [2, 2]
    });
    var modGridyAxis = yAxis.get("renderer");
    modGridyAxis.grid.template.setAll({
        strokeDasharray: [2, 2]
    });


    // Create series
    // https://www.amcharts.com/docs/v5/charts/xy-chart/series/
    var series = chart.series.push(am5xy.ColumnSeries.new(root, {
        xAxis: xAxis,
        yAxis: yAxis,
        valueYField: "value",
        categoryXField: "category",
        stroke: am5.color("#8E8A8A"),
        adjustBulletPosition: true,
        tooltip: am5.Tooltip.new(root, {
            labelText: "{valueY}"
        })
        }));

        series.columns.template.setAll({
        width: 0.5
        });

    //     series.bullets.push(function() {
    //         return am5.Bullet.new(root, {
    //             locationY: 1,
    //             sprite: am5.Circle.new(root, {
    //             radius: 5,
    //             fill: am5.color(0x000000),
                
    //             })
    //         })
    // })

    series.bullets.push(function(_root, _series, dataItem) {
        var bulletCircle = am5.Circle.new(root, {
          radius: 5,
          locationY: 1,
          fill: am5.color("#8E8A8A"),
          fillOpacity: 1,
        });
      
    bulletCircle.adapters.add("fill", function(fill, target) {
        var dataItem = target.dataItem;
        if (dataItem) {
        var sum = dataItem.dataContext.value
        var min = Infinity;
        var max = -Infinity;
        am5.array.each(series.dataItems, function(dataItem) {
            var value = dataItem.dataContext.value;
            if (value < min) {
            min = value;
            }
            if (value > max) {
            max = value;
            }
        })
    
        return am5.Color.interpolate((sum - min) / (max - min), am5.color(0x673AB7), am5.color(0xF44336));
        }
        return fill;
    })
    
    return am5.Bullet.new(root, {
        sprite: bulletCircle, 
        locationY: 1,
        // radius: 5
    });
    });
    
    series.children.moveValue(series.bulletsContainer, series.children.length + 1);


    series.on("tooltipDataItem", function (tooltipDataItem) {

        var residueID = tooltipDataItem.dataContext.category.split(' ')[1]
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


    series.set("heatRules", [{
        target: series.columns.template,
        key: "fill",
        min: am5.color(0x673AB7),
        max: am5.color(0xF44336),
        dataField: "valueY",
        key: "fill"
    }]);

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