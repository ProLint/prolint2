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
root.setThemes([am5themes_Animated.new(root), myTheme]);

// tell that valueX should be formatted as a date (show week number)
root.dateFormatter.setAll({
    dateFormat: "w",
    dateFields: ["valueX"]
});

root.locale.firstDayOfWeek = 0;

// var data = {};
fetch('/data/' + document.getElementById('lipids').value)
    .then(response => response.json())
    .then(data => {

        function prepareDataInput(data) {

            var weeklyData = [];
            var dailyData = [];
            var weekdays = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
            var weekAxisData = [
                { day: "Sun" },
                { day: "Mon" },
                { day: "Tue" },
                { day: "Wed" },
                { day: "Thu" },
                { day: "Fri" },
                { day: "Sat" }
                ];
       
            var firstDay = am5.time.round(new Date(data[0]["Activity Date"]), "year", 1); // Tue Jan 01 2019 00:00:00 GMT-0600 (Central Standard Time)
            // console.log('firstDay', firstDay);
            var total = 0;
            for (var i = 0; i < 53; i++) {
                weeklyData[i] = {};
                weeklyData[i].distance = 0;
                var date = new Date(firstDay);
                date.setDate(i * 7);
                am5.time.round(date, "week", 1);
                var endDate = am5.time.round(new Date(date), "week", 1);

                weeklyData[i].date = date.getTime();
                weeklyData[i].endDate = endDate.getTime();
            }

            am5.array.each(data, function (di) {
                var date = new Date(di["Activity Date"]);
                var weekDay = date.getDay();
                var weekNumber = am5.utils.getWeek(date);

                if (weekNumber == 1 && date.getMonth() == 11) {
                    weekNumber = 53;
                }

                var distance = am5.math.round(di["Distance"] / 1000, 1);

                weeklyData[weekNumber - 1].distance += distance;
                weeklyData[weekNumber - 1].distance = am5.math.round(
                    weeklyData[weekNumber - 1].distance,
                    1
                );
                total += distance;

                dailyData.push({
                    date: date.getTime(),
                    day: weekdays[weekDay],
                    "Distance": distance,
                    title: di["Activity Name"]
                });
            });

            return {
                dailyData,
                weeklyData,
                total,
                firstDay,
                weekAxisData
            };
        }

        var colorSet = am5.ColorSet.new(root, {});

        preparedData = prepareDataInput(data);

        var firstDay = preparedData['firstDay'];
        var weekAxisData = preparedData['weekAxisData'];

        console.log(preparedData);

        // Create chart
        // https://www.amcharts.com/docs/v5/charts/radar-chart/
        var chart = root.container.children.push(
            am5radar.RadarChart.new(root, {
                panX: false,
                panY: false,
                wheelX: "panX",
                wheelY: "zoomX",
                innerRadius: am5.percent(38),
                radius: am5.percent(85),
                startAngle: 270 - 170,
                endAngle: 270 + 170
                // startAngle: 270 - 170,
                // endAngle: 290 + 170
            })
        );

        // add label in the center
        chart.radarContainer.children.push(
            am5.Label.new(root, {
                text: "[fontSize:1.5em]ProLint[/]\n[fontSize:0.7em]" +
                    Math.random() +
                    " km[/]",
                textAlign: "center",
                centerX: am5.percent(50),
                centerY: am5.percent(50)
            })
        );

        // Add cursor
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Cursor
        var cursor = chart.set(
            "cursor",
            am5radar.RadarCursor.new(root, {
                behavior: "zoomX"
            })
        );
        cursor.lineY.set("visible", false);

        // Create axes and their renderers
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Adding_axes

        // date axis
        var dateAxisRenderer = am5radar.AxisRendererCircular.new(root, {
            minGridDistance: 20
        });

        dateAxisRenderer.labels.template.setAll({
            radius: 30,
            textType: "radial",
            centerY: am5.p50
        });

        // TODO: Need to change this to something ResName:ResNumber
        var dateAxis = chart.xAxes.push(
            am5xy.DateAxis.new(root, {
                baseInterval: {
                    timeUnit: "week",
                    count: 1
                },
                renderer: dateAxisRenderer,
                min: new Date(2019, 0, 1, 0, 0, 0).getTime(),
                max: new Date(2020, 0, 1, 0, 0, 0).getTime()
            })
        );

        // distance axis
        var distanceAxisRenderer = am5radar.AxisRendererRadial.new(root, {
            // axisAngle: 90,
            // radius: am5.percent(60),
            // innerRadius: am5.percent(20),
            axisAngle: 90,
            innerRadius: am5.percent(70),
            radius: am5.percent(100),

            // inversed: true,
            minGridDistance: 20
        });

        distanceAxisRenderer.labels.template.setAll({
            centerX: am5.p50,
            minPosition: 0.05,
            maxPosition: 0.95
        });

        var distanceAxis = chart.yAxes.push(
            am5xy.ValueAxis.new(root, {
                renderer: distanceAxisRenderer
            })
        );

        distanceAxis.set("numberFormat", "# ' km'");

        // week axis
        var weekAxisRenderer = am5radar.AxisRendererRadial.new(root, {
            axisAngle: 90,
            radius: am5.percent(68),
            innerRadius: am5.percent(40),
            // axisAngle: 90,
            // innerRadius: am5.percent(60),
            // radius: am5.percent(100),
            minGridDistance: 20
        });

        weekAxisRenderer.labels.template.setAll({
            centerX: am5.p50
        });

        var weekAxis = chart.yAxes.push(
            am5xy.CategoryAxis.new(root, {
                categoryField: "day",
                renderer: weekAxisRenderer
            })
        );

        // Create series
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Adding_series
        var distanceSeries = chart.series.push(
            am5radar.RadarColumnSeries.new(root, {
                calculateAggregates: true,
                xAxis: dateAxis,
                yAxis: distanceAxis,
                valueYField: "distance",
                valueXField: "date",
                tooltip: am5.Tooltip.new(root, {
                    labelText: "week {valueX}: {valueY}"
                })
            })
        );

        distanceSeries.columns.template.set("strokeOpacity", 0);

        // Set up heat rules
        // https://www.amcharts.com/docs/v5/concepts/settings/heat-rules/
        distanceSeries.set("heatRules", [{
            target: distanceSeries.columns.template,
            key: "fill",
            min: am5.color(0x673ab7),
            max: am5.color(0xf44336),
            dataField: "valueY"
        }]);

        // bubble series is a line series with stroeks hiddden
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Adding_series
        var bubbleSeries = chart.series.push(
            am5radar.RadarLineSeries.new(root, {
                calculateAggregates: true,
                xAxis: dateAxis,
                yAxis: weekAxis,
                baseAxis: dateAxis,
                categoryYField: "day",
                valueXField: "date",
                valueField: "Distance",
                maskBullets: false
            })
        );

        // only bullets are visible, hide stroke
        bubbleSeries.strokes.template.set("forceHidden", true);

        // add bullet
        var circleTemplate = am5.Template.new({});
        bubbleSeries.bullets.push(function () {
            var graphics = am5.Circle.new(root, {
                fill: distanceSeries.get("fill"),
                tooltipText: "{title}: {value} km"
            }, circleTemplate);
            return am5.Bullet.new(root, {
                sprite: graphics
            });
        });

        // Add heat rule (makes bubbles to be of a various size, depending on a value)
        // https://www.amcharts.com/docs/v5/concepts/settings/heat-rules/
        bubbleSeries.set("heatRules", [{
            target: circleTemplate,
            min: 3,
            max: 15,
            dataField: "value",
            key: "radius"
        }]);

        // set data
        // https://www.amcharts.com/docs/v5/charts/radar-chart/#Setting_data
        distanceSeries.data.setAll(preparedData['weeklyData']);
        weekAxis.data.setAll(weekAxisData);
        bubbleSeries.data.setAll(preparedData['dailyData']);

        bubbleSeries.appear(1000);
        distanceSeries.appear(1000);
        chart.appear(1000, 100);

        // create axis ranges
        var months = [
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec"
        ];
        for (var i = 0; i < 12; i++) {
            createRange(months[i], i);
        }

        function createRange(name, index) {
            var axisRange = dateAxis.createAxisRange(
                dateAxis.makeDataItem({
                    above: true
                })
            );
            axisRange.get("label").setAll({
                text: name
            });

            var fromTime = new Date(firstDay.getFullYear(), i, 1, 0, 0, 0).getTime();
            var toTime = am5.time.add(new Date(fromTime), "month", 1).getTime();

            axisRange.set("value", fromTime);
            axisRange.set("endValue", toTime);

            // every 2nd color for a bigger contrast
            var fill = axisRange.get("axisFill");
            fill.setAll({
                toggleKey: "active",
                cursorOverStyle: "pointer",
                fill: colorSet.getIndex(index * 2),
                visible: true,
                dRadius: 25,
                innerRadius: -25
            });
            axisRange.get("grid").set("visible", false);

            var label = axisRange.get("label");
            label.setAll({
                fill: am5.color(0xffffff),
                textType: "circular",
                radius: 8,
                text: months[index]
            });

            // clicking on a range zooms in
            fill.events.on("click", function (event) {
                var dataItem = event.target.dataItem;
                if (event.target.get("active")) {
                    dateAxis.zoom(0, 1);
                } else {
                    dateAxis.zoomToValues(dataItem.get("value"), dataItem.get("endValue"));
                }
            });
        }

        document.getElementById('lipids').addEventListener('change', function (e) {

            fetch('/data/' + e.target.value)
                .then(response => response.json())
                .then(data => {
                    preparedData = prepareDataInput(data);
                    distanceSeries.data.setAll(preparedData['weeklyData']);
                    weekAxis.data.setAll(weekAxisData);
                    bubbleSeries.data.setAll(preparedData['dailyData']);
                });
        });

    });