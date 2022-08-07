/**
 * ---------------------------------------
 * This app was created using amCharts 5.
 *
 * For more information visit:
 * https://www.amcharts.com/
 *
 * Documentation is available at:
 * https://www.amcharts.com/docs/v5/
 * ---------------------------------------
 */

var root = am5.Root.new("chartdiv");

const theme = am5.Theme.new(root);
theme.rule("Label").set("fontSize", 10);
theme.rule("Grid").set("strokeOpacity", 0.06);

root.setThemes([
    am5themes_Animated.new(root),
    theme
]);

root.fps = 60;
// Data
// fetch('/data/' + document.getElementById('lipids').value)
fetch('girk.json')
    .then(response => response.json())
    .then(contactData => {

        // console.log('contactData', contactData)

        var startFrameGroup = 0;
        var endFrameGroup = 1;
        var currentFrameGroup = 1;

        var div = document.getElementById("chartdiv");

        var colorSet = am5.ColorSet.new(root, {});

        // Params
        var innerRadius = 20;
        var lipidSelection = "CHOL";

        // Create chart
        // https://www.amcharts.com/docs/v5/charts/radar-chart/
        var chart = root.container.children.push(am5radar.RadarChart.new(root, {
            panX: false,
            panY: false,
            wheelX: "panX",
            wheelY: "zoomX",
            innerRadius: am5.percent(innerRadius),
            radius: am5.percent(65),
            startAngle: 270 - 170,
            endAngle: 270 + 170
            // Right handed half circle:
            // startAngle: -90,
            // endAngle: 90
        }));


        // Add cursor
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Cursor
        var cursor = chart.set("cursor", am5radar.RadarCursor.new(root, {
            behavior: "zoomX",
            radius: am5.percent(innerRadius),
            innerRadius: -25
        }));
        cursor.lineY.set("visible", true);
        cursor.lineY.set("opacity", 0.5);

        // Create axes and their renderers
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Adding_axes
        var xRenderer = am5radar.AxisRendererCircular.new(root, {
            minGridDistance: 10
        });

        xRenderer.labels.template.setAll({
            radius: 10,
            textType: "radial",
            centerY: am5.p50
        });

        var yRenderer = am5radar.AxisRendererRadial.new(root, {
            axisAngle: 90
        });

        yRenderer.labels.template.setAll({
            centerX: am5.p50
        });

        var categoryAxis = chart.xAxes.push(am5xy.CategoryAxis.new(root, {
            maxDeviation: 0,
            categoryField: "residue",
            renderer: xRenderer
        }));

        var valueAxis = chart.yAxes.push(am5xy.ValueAxis.new(root, {
            min: 0,
            max: 4,
            extraMax: 0.1,
            renderer: yRenderer
        }));


        // Create series
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Adding_series
        var series = chart.series.push(am5radar.RadarColumnSeries.new(root, {
            calculateAggregates: true,
            name: "Series",
            xAxis: categoryAxis,
            yAxis: valueAxis,
            valueYField: "value_" + currentFrameGroup,
            categoryXField: "residue",
            tooltip: am5.Tooltip.new(root, {
                labelText: "{categoryX}: {valueY}"
            })
        }));

        series.columns.template.set("strokeOpacity", 0);


        // Set up heat rules
        // https://www.amcharts.com/docs/v5/concepts/settings/heat-rules/
        series.set("heatRules", [{
            target: series.columns.template,
            key: "fill",
            min: am5.color(0x673AB7),
            max: am5.color(0xF44336),
            dataField: "valueY"
        }]);

        // Add scrollbars
        // https://www.amcharts.com/docs/v5/charts/xy-chart/scrollbars/
        // chart.set("scrollbarX", am5.Scrollbar.new(root, {
        //     orientation: "horizontal"
        // }));
        // chart.set("scrollbarY", am5.Scrollbar.new(root, {
        //     orientation: "vertical"
        // }));

        // Add year label
        var yearLabel = chart.radarContainer.children.push(am5.Label.new(root, {
            fontSize: "2em",
            text: currentFrameGroup.toString(),
            centerX: am5.p50,
            centerY: am5.p50,
            fill: am5.color(0x673AB7)
        }));


        // Generate and set data
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Setting_data
        var data = generateRadarData();
        series.data.setAll(data);
        categoryAxis.data.setAll(data);

        series.appear(10);
        chart.appear(10, 100);

        function generateRadarData() {
            var data = [];
            var i = 0;
            for (var lipid in contactData) {
                if (lipid != lipidSelection) {
                    continue
                }
                var lipidData = contactData[lipid];

                lipidData.forEach(function (residue) {
                    var rawDataItem = {
                        "residue": residue[0]
                    }

                    var startY = 1
                    for (var y = startY; y < residue.length; y++) {
                        rawDataItem["value_" + (startFrameGroup + y - startY)] = residue[y];
                    }

                    data.push(rawDataItem);
                });

                createRange(lipid, lipidData, i);
                i++;

            }
            // console.log('dataProcessed', data)
            return data;
        }


        function createRange(name, lipidData, index) {
            var axisRange = categoryAxis.createAxisRange(categoryAxis.makeDataItem({
                above: true
            }));
            axisRange.get("label").setAll({
                text: name
            });
            // first residue
            axisRange.set("category", lipidData[0][0]);
            // last residue
            axisRange.set("endCategory", lipidData[lipidData.length - 1][0]);

            // every 3rd color for a bigger contrast
            var fill = axisRange.get("axisFill");
            fill.setAll({
                toggleKey: "active",
                cursorOverStyle: "pointer",
                fill: colorSet.getIndex(index * 3),
                // fill: colorSet.getIndex(3),
                visible: true,
                innerRadius: -25
            });
            axisRange.get("grid").set("visible", false);

            var label = axisRange.get("label");
            label.setAll({
                fill: am5.color(0xffffff),
                textType: "circular",
                visible: true,
                radius: -16
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


        // Create controls
        var container = chart.children.push(am5.Container.new(root, {
            y: am5.percent(95),
            centerX: am5.p50,
            x: am5.p50,
            width: am5.percent(40),
            layout: root.horizontalLayout
        }));

        var playButton = container.children.push(am5.Button.new(root, {
            themeTags: ["play"],
            visible: false,
            centerY: am5.p50,
            marginRight: 15,
            icon: am5.Graphics.new(root, {
                themeTags: ["icon"]
            })
        }));

        playButton.events.on("click", function () {
            if (playButton.get("active")) {
                slider.set("start", slider.get("start") + 0.0001);
            } else {
                slider.animate({
                    key: "start",
                    to: 1,
                    duration: 15000 * (1 - slider.get("start"))
                });
            }
        })

        var slider = container.children.push(am5.Slider.new(root, {
            orientation: "horizontal",
            visible: false,
            start: 0.0,
            centerY: am5.p50
        }));

        slider.on("start", function (start) {
            if (start === 1) {
                playButton.set("active", false);
            }
        });

        slider.events.on("rangechanged", function () {
            // val = Math.round(slider.get("start", 0) * (endFrameGroup - startFrameGroup));
            // val = slider.get("start", 0) //* (endFrameGroup - startFrameGroup)
            // console.log('before UPDATE', val)
            updateRadarData(startFrameGroup + Math.round(slider.get("start", 0) * (endFrameGroup - startFrameGroup)));
        });

        function updateRadarData(year) {
            if (currentFrameGroup != year) {
                // console.log('INSIDE UPDATE', year)
                currentFrameGroup = year;
                yearLabel.set("text", currentFrameGroup.toString());
                am5.array.each(series.dataItems, function (dataItem) {
                    var newValue = dataItem.dataContext["value_" + year];
                    dataItem.set("valueY", newValue);
                    dataItem.animate({
                        key: "valueYWorking",
                        to: newValue,
                        duration: 500
                    });
                });
            }
        }

    });