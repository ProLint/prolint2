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

// Create root element
// https://www.amcharts.com/docs/v5/getting-started/#Root_element
var root = am5.Root.new("chartdiv");

// Create custom theme
// https://www.amcharts.com/docs/v5/concepts/themes/#Quick_custom_theme
const myTheme = am5.Theme.new(root);
myTheme.rule("Label").set("fontSize", 10);
myTheme.rule("Grid").set("strokeOpacity", 0.06);

// Set themes
// https://www.amcharts.com/docs/v5/concepts/themes/
root.setThemes([
    am5themes_Animated.new(root),
    myTheme
]);

// Data
// fetch('/data/' + document.getElementById('lipids').value)
fetch('d.json')
    .then(response => response.json())
    .then(temperatures => {

        console.log('temps', temperatures)

        // Modify defaults
        root.numberFormatter.set("numberFormat", "+#.0°C|#.0°C|0.0°C");

        var startYear = 1973;
        var endYear = 2016;
        var currentYear = 1974;

        var div = document.getElementById("chartdiv");

        var colorSet = am5.ColorSet.new(root, {});

        // Params
        var innerRadius = 20

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
            categoryField: "country",
            renderer: xRenderer
        }));

        var valueAxis = chart.yAxes.push(am5xy.ValueAxis.new(root, {
            min: -3,
            max: 6,
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
            valueYField: "value" + currentYear,
            categoryXField: "country",
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
            text: currentYear.toString(),
            centerX: am5.p50,
            centerY: am5.p50,
            fill: am5.color(0x673AB7)
        }));


        // Generate and set data
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Setting_data
        var data = generateRadarData();
        console.log('data', data)
        series.data.setAll(data);
        categoryAxis.data.setAll(data);

        series.appear(500);
        chart.appear(500, 100);

        function generateRadarData() {
            var data = [];
            var i = 0;
            for (var continent in temperatures) {
                var continentData = temperatures[continent];

                continentData.forEach(function (country) {
                    var rawDataItem = {
                        "country": country[0]
                    }

                    for (var y = 2; y < country.length; y++) {
                        rawDataItem["value" + (startYear + y - 2)] = country[y];
                    }

                    data.push(rawDataItem);
                });

                createRange(continent, continentData, i);
                i++;

            }
            return data;
        }


        function createRange(name, continentData, index) {
            var axisRange = categoryAxis.createAxisRange(categoryAxis.makeDataItem({
                above: true
            }));
            axisRange.get("label").setAll({
                text: name
            });
            // first country
            axisRange.set("category", continentData[0][0]);
            // last country
            axisRange.set("endCategory", continentData[continentData.length - 1][0]);

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
            start: 0.0,
            centerY: am5.p50
        }));

        slider.on("start", function (start) {
            if (start === 1) {
                playButton.set("active", false);
            }
        });

        slider.events.on("rangechanged", function () {
            // val = Math.round(slider.get("start", 0)* (endYear - startYear));
            // val = slider.get("start", 0) //* (endYear - startYear)
            // console.log('before UPDATE', val)
            updateRadarData(startYear + Math.round(slider.get("start", 0) * (endYear - startYear)));
        });

        function updateRadarData(year) {
            if (currentYear != year) {
                // console.log('INSIDE UPDATE', year)
                currentYear = year;
                yearLabel.set("text", currentYear.toString());
                am5.array.each(series.dataItems, function (dataItem) {
                    var newValue = dataItem.dataContext["value" + year];
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