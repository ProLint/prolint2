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


// Fetch the data from the backend
var obj = {
    "lipid": "",
    "protein": ""
}
fetch('/data/' + JSON.stringify(obj))
    .then(response => response.json())
    .then(responseData => {

        // console.log('responseData', responseData);
        // console.log('top_lipids', responseData['globalTopLipids'])
        // console.log('lipid_contact_frames', responseData['lipidContactFrames'])
        var contactData = responseData['data'];
        var lipids = responseData['lipids'];
        var proteins = responseData['proteins'];

        var systemHasOneProtein = false;
        if (proteins.length == 1) {
            systemHasOneProtein = true
        }

        ///////////////////////////////////////////
        /////////////// RadarApp //////////////////
        ///////////////////////////////////////////
        var root = am5.Root.new("chartdiv");

        const theme = am5.Theme.new(root);
        theme.rule("Label").set("fontSize", 10);
        theme.rule("Grid").set("strokeOpacity", 0.06);

        root.setThemes([
            am5themes_Animated.new(root),
            theme
        ]);

        root.fps = 60;

        var startFrameGroup = 0;
        var endFrameGroup = 1;
        var currentFrameGroup = 1;

        var colorSet = am5.ColorSet.new(root, {});

        // Params
        var innerRadius = 50;

        // Create chart
        var chart = root.container.children.push(am5radar.RadarChart.new(root, {
            panX: false,
            panY: false,
            wheelX: "panX",
            wheelY: "zoomX",
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
            innerRadius: -25
        }));
        cursor.lineY.set("visible", true);
        cursor.lineY.set("opacity", 0.5);

        // Create axes and their renderers
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
            strictMinMax: true,
            extraMax: 0.1,
            renderer: yRenderer
        }));

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
            valueYField: "value_" + currentFrameGroup,
            categoryXField: "residue",
            tooltip: am5.Tooltip.new(root, {
                labelText: "{categoryX}: {valueY}"
            })
        }));

        series.columns.template.set("strokeOpacity", 0);

        // Set up heat rules
        series.set("heatRules", [{
            target: series.columns.template,
            key: "fill",
            min: am5.color(0x673AB7),
            max: am5.color(0xF44336),
            dataField: "valueY"
        }]);

        // Add scrollbars
        // chart.set("scrollbarX", am5.Scrollbar.new(root, {
        //     orientation: "horizontal"
        // }));
        // chart.set("scrollbarY", am5.Scrollbar.new(root, {
        //     orientation: "vertical"
        // }));

        // Add frameGroup label
        // var frameGroupLabel = chart.radarContainer.children.push(am5.Label.new(root, {
        //     fontSize: "2em",
        //     text: currentFrameGroup.toString(),
        //     centerX: am5.p50,
        //     centerY: am5.p50,
        //     fill: am5.color(0x673AB7)
        // }));

        // Generate and set data
        var data = generateRadarData(contactData);
        series.data.setAll(data);
        categoryAxis.data.setAll(data);

        series.appear(500);
        chart.appear(500, 100);

        function generateRadarData(contactData) {
            // contactData = contactData['Protein0']
            var data = [];
            var i = 0;
            for (var lipid in contactData) {
                var lipidData = contactData[lipid];

                lipidData.forEach(function (residue) {
                    var rawDataItem = {
                        "residue": residue[0]
                    }

                    var startY = 1
                    for (var y = startY; y < residue.length; y++) {
                        rawDataItem["value_" + (startFrameGroup + y - startY)] = residue[y];
                    }
                    // rawDataItem['protein'] = "GIRK"
                    data.push(rawDataItem);
                });

                createRange(lipid, lipidData, i);
                i++;

            }
            return data;
        }

        function createRange(name, lipidData, index) {
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
            // console.log('before UPDATE')
            updateRadarData(startFrameGroup + Math.round(slider.get("start", 0) * (endFrameGroup - startFrameGroup)));
        });

        function updateRadarData(frameGroup) {
            if (currentFrameGroup != frameGroup) {
                currentFrameGroup = frameGroup;
                // frameGroupLabel.set("text", currentFrameGroup.toString());
                am5.array.each(series.dataItems, function (dataItem) {
                    var newValue = dataItem.dataContext["value_" + frameGroup];
                    dataItem.set("valueY", newValue);
                    dataItem.animate({
                        key: "valueYWorking",
                        to: newValue,
                        duration: 500
                    });
                });
            }
        }

        ///////////////////////////////////////////
        /////////////// PieApp ////////////////////
        ///////////////////////////////////////////
        var pieRoot = am5.Root.new("chartdiv2");

        // Set themes
        pieRoot.setThemes([am5themes_Animated.new(pieRoot)]);

        var pieContainer = pieRoot.container.children.push(
            am5.Container.new(pieRoot, {
                width: am5.p100,
                height: am5.p100,
                layout: pieRoot.horizontalLayout
            })
        );

        // Create main chart
        var pieChart = pieContainer.children.push(
            am5percent.PieChart.new(pieRoot, {
                tooltip: am5.Tooltip.new(pieRoot, {})
            })
        );

        // Create series
        var pieSeries = pieChart.series.push(
            am5percent.PieSeries.new(pieRoot, {
                valueField: "value",
                categoryField: "category",
                alignLabels: false
            })
        );

        pieSeries.labels.template.setAll({
            textType: "circular",
            radius: 4
        });
        pieSeries.ticks.template.set("visible", false);
        pieSeries.slices.template.set("toggleKey", "none");

        // add events
        pieSeries.slices.template.events.on("click", function (e) {
            selectSlice(e.target);

            // console.log('series', series)

            // TODO:
            // only execute when the protein changes
            var protein = e.target.dataItem.dataContext.category
            var lipid = subSeries.slices.getIndex(0).dataItem.dataContext.category

            obj.lipid = lipid
            obj.protein = protein
            fetch('/data/' + JSON.stringify(obj))
                .then(response => response.json())
                .then(responseData => {

                    updateData = responseData['data']
                    var updateData = generateRadarData(updateData);
                    series.data.setAll(updateData);
                    // categoryAxis.data.setAll(updateData);

                    am5.array.each(series.dataItems, function (dataItem) {
                        var newValue = dataItem.dataContext["value_" + 0];
                        dataItem.set("valueY", newValue);
                        dataItem.animate({
                            key: "valueYWorking",
                            to: newValue,
                            duration: 0
                        });
                    });
                });
            series.appear(1000);
            // chart.appear(500, 100);

        });

        // Create sub chart
        var subChart = pieContainer.children.push(
            am5percent.PieChart.new(pieRoot, {
                radius: am5.percent(50),
                tooltip: am5.Tooltip.new(pieRoot, {})
            })
        );

        // Create sub series
        var subSeries = subChart.series.push(
            am5percent.PieSeries.new(pieRoot, {
                valueField: "value",
                categoryField: "category",
                alignLabels: false
            })
        );

        // subSeries click event to link to radar chart
        subSeries.slices.template.events.on("click", function (e) {
            var lipid = e.target.dataItem.dataContext.category
            if (lipid != axisRange.get("label").get('text')) {
                obj.lipid = lipid
                // TODO:
                // get correct protein
                obj.protein = "GIRK"
                fetch('/data/' + JSON.stringify(obj))
                    .then(response => response.json())
                    .then(responseData => {

                        updateData = responseData['data']

                        var updateData = generateRadarData(updateData);
                        series.data.setAll(updateData);
                        // categoryAxis.data.setAll(updateData);

                        am5.array.each(series.dataItems, function (dataItem) {
                            var newValue = dataItem.dataContext["value_" + 0];
                            dataItem.set("valueY", newValue);
                            dataItem.animate({
                                key: "valueYWorking",
                                to: newValue,
                                duration: 0
                            });
                        });
                    });
                series.appear(1000);
                // chart.appear(500, 100);
            }
        });

        subSeries.data.setAll(lipids.map(lipidName => ({
            category: lipidName,
            value: 0
        })))
        subSeries.labels.template.setAll({
            textType: "circular",
            radius: 4
        });
        subSeries.ticks.template.set("visible", false);
        subSeries.slices.template.set("toggleKey", "none");

        var selectedSlice;

        pieSeries.on("startAngle", function () {
            updateLines();
        });

        //   pieContainer.events.on("boundschanged", function() {
        //     pieRoot.events.on("frameended", function(){
        //       updateLines();
        //      })
        //   })

        function updateLines() {
            if (selectedSlice) {
                var startAngle = selectedSlice.get("startAngle");
                var arc = selectedSlice.get("arc");
                var radius = selectedSlice.get("radius");

                if (!systemHasOneProtein) {
                    var x00 = radius * am5.math.cos(startAngle);
                    var y00 = radius * am5.math.sin(startAngle);

                    var x10 = radius * am5.math.cos(startAngle + arc);
                    var y10 = radius * am5.math.sin(startAngle + arc);

                } else {
                    var x00 = radius * am5.math.sin(startAngle);
                    var y00 = radius * am5.math.cos(startAngle);

                    var x10 = radius * am5.math.sin(startAngle + arc);
                    var y10 = -radius * am5.math.cos(startAngle + arc);

                }

                var subRadius = subSeries.slices.getIndex(0).get("radius");
                var x01 = 0;
                var y01 = -subRadius;

                var x11 = 0;
                var y11 = subRadius;

                var point00 = pieSeries.toGlobal({
                    x: x00,
                    y: y00
                });
                var point10 = pieSeries.toGlobal({
                    x: x10,
                    y: y10
                });

                var point01 = subSeries.toGlobal({
                    x: x01,
                    y: y01
                });
                var point11 = subSeries.toGlobal({
                    x: x11,
                    y: y11
                });

                line0.set("points", [point00, point01]);
                line1.set("points", [point10, point11]);
            }
        }

        // lines
        var line0 = pieContainer.children.push(
            am5.Line.new(pieRoot, {
                position: "absolute",
                stroke: pieRoot.interfaceColors.get("text"),
                strokeDasharray: [2, 2]
            })
        );
        var line1 = pieContainer.children.push(
            am5.Line.new(pieRoot, {
                position: "absolute",
                stroke: pieRoot.interfaceColors.get("text"),
                strokeDasharray: [2, 2]
            })
        );

        // Set data
        pieSeries.data.setAll(responseData['pieData']);

        function selectSlice(slice) {
            selectedSlice = slice;
            var dataItem = slice.dataItem;
            var dataContext = dataItem.dataContext;

            if (dataContext) {
                var i = 0;
                subSeries.data.each(function (dataObject) {
                    var dataObj = dataContext.subData[i];
                    if (dataObj) {
                        subSeries.data.setIndex(i, dataObj);
                        if (!subSeries.dataItems[i].get("visible")) {
                            subSeries.dataItems[i].show();
                        }
                    } else {
                        subSeries.dataItems[i].hide();
                    }

                    i++;
                });
            }

            var middleAngle = slice.get("startAngle") + slice.get("arc") / 2;
            var firstAngle = pieSeries.dataItems[0].get("slice").get("startAngle");

            pieSeries.animate({
                key: "startAngle",
                to: firstAngle - middleAngle,
                duration: 1000,
                easing: am5.ease.out(am5.ease.cubic)
            });
            pieSeries.animate({
                key: "endAngle",
                to: firstAngle - middleAngle + 360,
                duration: 1000,
                easing: am5.ease.out(am5.ease.cubic)
            });
        }

        pieContainer.appear(1000, 10);

        pieSeries.events.on("datavalidated", function () {
            selectSlice(pieSeries.slices.getIndex(0));
        });



        ///////////////////////////////////////////
        ////////////// Lipid Table ////////////////
        ///////////////////////////////////////////
        var table = new Tabulator("#lipid-table", {
            data: responseData['tableData'],
            height: "300px",
            columns: [{
                    title: "Lipid ID",
                    field: "lipidID",
                    width: 70,
                    hozAlign:"center"
                },
                {
                    title: "Total Contacts",
                    field: "contactFrequency",
                    width: 70,
                    hozAlign:"center"
                },
            ],
            headerVisible: false,
        });

        table.on("rowClick", function (e, row) {

            obj = {
                "lipidID": row.getData()['lipidID']
            }
            fetch('/toplipids/' + JSON.stringify(obj))
                .then(response => response.json())
                .then(tableResponseData => {

                    ganttData = tableResponseData['ganttData'].map((lp, ix) => ({
                        ...lp,
                        columnSettings: {
                            fill: colorSet.getIndex(ix * 3)
                        }
                    }))
                    ganttYAxis.data.setAll(tableResponseData['topLipids'].map(v => ({
                        category: v
                    })))
                    ganttSeries.data.setAll(ganttData);


                });

        });

        ///////////////////////////////////////////
        /////////////// GanttApp //////////////////
        ///////////////////////////////////////////
        var ganttRoot = am5.Root.new("chartdiv3");

        // Set themes
        ganttRoot.setThemes([
            am5themes_Animated.new(ganttRoot)
        ]);


        // Create chart
        var ganttChart = ganttRoot.container.children.push(am5xy.XYChart.new(ganttRoot, {
            panX: false,
            panY: false,
            wheelX: "panX",
            wheelY: "zoomX",
            layout: ganttRoot.verticalLayout
        }));

        var legend = ganttChart.children.push(am5.Legend.new(ganttRoot, {
            centerX: am5.p50,
            x: am5.p50
        }))

        var colors = ganttChart.get("colors");

        // Data
        ganttData = responseData['ganttData'].map((lp, ix) => ({
            ...lp,
            columnSettings: {
                fill: colorSet.getIndex(ix * 3)
            }
        }))

        // Create axes
        var ganttYAxis = ganttChart.yAxes.push(
            am5xy.CategoryAxis.new(ganttRoot, {
                categoryField: "category",
                renderer: am5xy.AxisRendererY.new(ganttRoot, {
                    inversed: true,
                    minGridDistance: 1
                }),
                tooltip: am5.Tooltip.new(ganttRoot, {
                    themeTags: ["axis"],
                    animationDuration: 200
                })
            })
        );

        // let ganttYRenderer = ganttYAxis.get("renderer");
        // ganttYRenderer.labels.template.setAll({
        //   fill: am5.color(0xFF0000),
        //   fontSize: "0.4em",
        // });
        ganttYAxis.data.setAll(responseData['topLipids'].map(v => ({
            category: v
        })))

        var ganttXAxis = ganttChart.xAxes.push(am5xy.ValueAxis.new(ganttRoot, {
            min: 0,
            max: 180,
            strictMinMax: true,
            renderer: am5xy.AxisRendererX.new(ganttRoot, {})
        }));

        var ganttSeries = ganttChart.series.push(am5xy.ColumnSeries.new(ganttRoot, {
            xAxis: ganttXAxis,
            yAxis: ganttYAxis,
            openValueXField: "startFrame",
            valueXField: "endFrame",
            categoryYField: "category",
            sequencedInterpolation: true
        }));

        ganttSeries.columns.template.setAll({
            templateField: "columnSettings",
            strokeOpacity: 0,
            fillOpacity: 0.8,
            tooltipText: "{category}"
        });

        ganttSeries.data.setAll(ganttData);

        ganttSeries.appear();
        ganttChart.appear(1000, 100);

        ganttSeries.columns.template.events.on("click", function (e, d) {
            residueID = e.target.dataItem.dataContext.category;

            ctx = e.target.dataItem.dataContext;

            obj = {
                "lipidID": ctx.lipid_id,
                "residueID": ctx.category
            }
            fetch('/distance/' + JSON.stringify(obj))
                .then(response => response.json())
                .then(heatmapResponseData => {
                    heatmapSeries.data.setAll(heatmapResponseData['heatmapData']);
                    hmYAxis.data.setAll(heatmapResponseData['residueAtomsData']);
                    hmXAxis.data.setAll(heatmapResponseData['lipidAtomsData']);
                });

        });

        ///////////////////////////////////////////
        ///////////// Heatmap App /////////////////
        ///////////////////////////////////////////
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
            min: am5.color(0xfffb77),
            max: am5.color(0xfe131a),
            dataField: "value",
            key: "fill"
        }]);

        // Add heat legend
        var heatLegend = heatmapChart.bottomAxesContainer.children.push(am5.HeatLegend.new(heatmapRoot, {
            orientation: "horizontal",
            endColor: am5.color(0xfffb77),
            startColor: am5.color(0xfe131a)
        }));

        // Set data
        heatmapSeries.data.setAll(responseData['heatmapData']);
        hmYAxis.data.setAll(responseData['residueAtomsData']);
        hmXAxis.data.setAll(responseData['lipidAtomsData']);

        // Make stuff animate on load
        heatmapChart.appear(1000, 100);

    });