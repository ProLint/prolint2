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

// NOTES (and TODO:)
// 1.
// The network app is always redrawn (rather than updated)
// This is mainly due to the pointerover events which retrieve the color
// value from the subseries pie chart. A solution might be update these
// events (or remove and add new events everytime the chart is updated).
// 2.
// The radarApp & networkApp use the same html element, but are rendered differently.
// The current solution is to use global objects to easily keep references to the disposed
// roots and releted elements. Ideally, we would not need to dispose root, but
// clear and reuse, and working within the same scope, to avoid using globals.

// 3. Data fetching is not optimized to the scope of the different objects. The JSON.stringify
// calls can also be avoided.

var rootReferenceObjects;
var networkRootReference;
var globalLipidReference;

networkRootReference = {
    "active": false
}

// Fetch the data from the backend
var obj = {
    "lipid": "",
    "protein": ""
}
fetch('/data/' + JSON.stringify(obj))
    .then(response => response.json())
    .then(responseData => {

        console.log('responseData', responseData)

        var frameNumber = responseData["frameNumber"];

        var contactData = responseData['data'];
        var lipids = responseData['lipids'];
        var proteins = responseData['proteins'];

        var startLipidID = responseData['tableData'][0].lipidID;
        var startResidueID = responseData['topLipids'][0]

        var systemHasOneProtein = false;
        if (proteins.length == 1) {
            systemHasOneProtein = true
        }

        ///////////////////////////////////////////
        /////////////// RadarApp //////////////////
        ///////////////////////////////////////////
        function radarApp() {

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
            viewerResidueSelection = {
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

                residueID = tooltipDataItem.dataContext.residue.split(' ')[1]
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

            // series.data.setAll(contactData);
            // categoryAxis.data.setAll(contactData);
            // createRange(lipid, contactData, 0);

            // series.appear(100);
            // chart.appear(100);

            var button = chart.plotContainer.children.push(am5.Button.new(root, {
                dx: 0,
                dy: 0,
                label: am5.Label.new(root, {
                    text: "Shared Contacts",
                    fontSize: 10,
                })
            }));

            button.events.on("click", function (ev) {
                root.dispose();

                am5.array.each(subSeries.dataItems, function (dataItem, ix) {
                    if (dataItem._settings.slice._settings.active) {
                        lipid = dataItem.dataContext.category
                        col = subSeries.get("colors").getIndex(ix)
                        networkRootReference = networkApp(lipid = lipid)
                        rootReferenceObjects["active"] = false
                    }
                })

            });

            return {
                chart: chart,
                series: series,
                axisRange: axisRange,
                categoryAxis: categoryAxis,
                createRange: createRange,
                active: true
            }
        }
        rootReferenceObjects = radarApp()
        rootReferenceObjects["series"].data.setAll(contactData);
        rootReferenceObjects["categoryAxis"].data.setAll(contactData);
        rootReferenceObjects["createRange"](lipids[0], contactData, 0);
        globalLipidReference = lipids[0]

        rootReferenceObjects["series"].appear(100);
        rootReferenceObjects["chart"].appear(100);

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
                innerRadius: am5.percent(80),
                tooltip: am5.Tooltip.new(pieRoot, {})
            })
        );

        // Create rootReferenceObjects["series"]
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

        var cols = pieSeries.get("colors")
        pieSeries.get("colors").set("colors", [
            am5.color(0x095256),
            am5.color(0x087f8c),
            am5.color(0x5aaa95),
            am5.color(0x86a873),
            am5.color(0xbb9f06)
        ]);

        pieSeries.labels.template.setAll({
            text: "{category}"
        });

        // add events
        pieSeries.slices.template.events.on("click", function (e) {
            selectSlice(e.target);

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
                    rootReferenceObjects["series"].data.setAll(updateData);
                    rootReferenceObjects["categoryAxis"].data.setAll(updateData);
                    rootReferenceObjects["createRange"](lipid, updateData, 0);

                });

            rootReferenceObjects["series"].appear(100);
        });

        // Create sub chart
        var subChart = pieContainer.children.push(
            am5percent.PieChart.new(pieRoot, {
                radius: am5.percent(50),
                innerRadius: am5.percent(50),
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
        var sliceTemplate = subSeries.slices.template;
        sliceTemplate.setAll({
            draggable: true,
            cornerRadius: 5,
            cursorOverStyle: "pointer",
        });

        sliceTemplate.events.on("pointerup", function (e) {
            var slice = e.target;

            subSeries.hideTooltip();
            slice.animate({
                key: "x",
                to: 0,
                duration: 500,
                easing: am5.ease.out(am5.ease.cubic)
            });
            slice.animate({
                key: "y",
                to: 0,
                duration: 500,
                easing: am5.ease.out(am5.ease.cubic)
            });
        });

        var currentSlice;
        subSeries.slices.template.on("active", function (active, slice) {
            if (currentSlice && currentSlice != slice && active) {
                currentSlice.set("active", undefined)
            }
            currentSlice = slice;
        });

        // subSeries click event to link to radar chart
        subSeries.slices.template.events.on("click", function (e) {

            // NOTE:
            // Along with the protein pieChart (not yet fully implemented),
            // this is the main entry point to the entire dashboard. Once a slice
            // is clicked, we need to update all apps of the dashboard. Below
            // we update the circular, table, gantt, and heatmap apps.
            var lipid_id = undefined;
            var residue_id = undefined;
            var lipid = e.target.dataItem.dataContext.category

            var sameLipid = true;
            if (rootReferenceObjects['active']) {
                if (lipid != rootReferenceObjects["axisRange"].get("label").get('text')) {
                    sameLipid = false;
                    // TODO:
                    // get correct protein
                    // Update Circular App Data
                    obj.protein = "GIRK"
                    obj.lipid = lipid
                    fetch('/data/' + JSON.stringify(obj))
                        .then(response => response.json())
                        .then(responseData => {

                            updateData = responseData['data'];

                            rootReferenceObjects["series"].data.setAll(updateData);
                            rootReferenceObjects["categoryAxis"].data.setAll(updateData);
                            rootReferenceObjects["createRange"](lipid, updateData, 0);

                            am5.array.each(subSeries.dataItems, function (dataItem, ix) {
                                if (dataItem.dataContext.category == lipid) {
                                    col = subSeries.get("colors").getIndex(ix)
                                    rootReferenceObjects["axisRange"].get("axisFill").set("fill", col)
                                }
                            })
                        });
                }
            } else {
                if (lipid != networkRootReference["series"]._settings.lipid) {
                    sameLipid = false;
                    networkRootReference["root"].dispose();
                    networkRootReference = networkApp(lipid = lipid)
                    networkRootReference["series"].appear(1000)

                }
            }

            if (!sameLipid) {

                // Update Table Data
                fetch('/tabledata/' + JSON.stringify(obj))
                    .then(response => response.json())
                    .then(pieChartResponseData => {
                        table.replaceData(pieChartResponseData['tableData']);
                        lipid_id = pieChartResponseData['tableData'][0]['lipidID']

                        // Update Gantt App Data
                        obj = {
                            "lipidID": lipid_id
                        }
                        fetch('/toplipids/' + JSON.stringify(obj))
                            .then(response => response.json())
                            .then(tableResponseData => {

                                ganttData = tableResponseData['ganttData'].map((lp, ix) => ({
                                    ...lp
                                }))
                                ganttYAxis.data.setAll(tableResponseData['topLipids'].map(v => ({
                                    category: v
                                })))
                                ganttSeries.data.setAll(ganttData);
                                sortCategoryAxis()

                                // Update Heatmap App Data
                                obj = {
                                    "lipidID": lipid_id,
                                    "residueID": ganttData[0]['category']
                                }
                                fetch('/distance/' + JSON.stringify(obj))
                                    .then(response => response.json())
                                    .then(heatmapResponseData => {
                                        heatmapSeries.data.setAll(heatmapResponseData['heatmapData']);
                                        hmYAxis.data.setAll(heatmapResponseData['residueAtomsData']);
                                        hmXAxis.data.setAll(heatmapResponseData['lipidAtomsData']);
                                    });

                            });

                    });

                rootReferenceObjects["series"].appear(100);
            }
        });

        subSeries.get("colors").set("colors", [
            am5.color(0x095256),
            am5.color(0x087f8c),
            am5.color(0x5aaa95),
            am5.color(0x86a873),
            am5.color(0xbb9f06)
        ]);

        subSeries.labels.template.setAll({
            text: "{category}"
        });

        subSeries.data.setAll(lipids.map(lipidName => ({
            category: lipidName,
            value: 0
        })))

        // subSeries.labels.template.setAll({
        //     textType: "circular",
        //     radius: 4
        // });
        // subSeries.ticks.template.set("visible", false);
        // subSeries.slices.template.set("toggleKey", "none");

        var selectedSlice;

        pieSeries.on("startAngle", function () {
            updateLines();
        });

        // Pre-select first slice
        subSeries.events.on("datavalidated", function () {
            subSeries.slices.getIndex(0).set("active", true);
        });
        pieContainer.events.on("boundschanged", function () {
            pieRoot.events.on("frameended", function () {
                updateLines();
            })
        })

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
            layout: "fitColumns",
            // autoResize:false,
            resizableRows: false,
            selectable: 1,
            selectablePersistence: false,
            columns: [{
                title: "Lipid <br>Contact Frequencies",
                columns: [{
                        title: "id",
                        field: "lipidID",
                        // width: 100,
                        hozAlign: "center",
                        // frozen:true,
                        headerSort: false,
                        resizable: false,
                        // headerFilter:"input"
                    },
                    {
                        title: "f",
                        field: "contactFrequency",
                        // width: 120,
                        hozAlign: "center",
                        headerSort: false,
                        resizable: false,
                    },
                ]
            }],
            headerVisible: true,
        });

        table.on("rowClick", function (e, row) {
            lipidID = row.getData()['lipidID']
            obj = {
                "lipidID": lipidID
            }
            fetch('/toplipids/' + JSON.stringify(obj))
                .then(response => response.json())
                .then(tableResponseData => {
                    ganttData = tableResponseData['ganttData'].map((lp, ix) => ({
                        ...lp,
                        besi: 'green'
                    }))
                    ganttYAxis.data.setAll(tableResponseData['topLipids'].map(v => ({
                        category: v
                    })))
                    ganttSeries.data.setAll(ganttData);
                    sortCategoryAxis()

                    // On LipidID selection, show interacting residue on the 3D viewer.
                    var selectSections = []
                    for (let ix = 0; ix < ganttData.length; ix++) {
                        const residueID = ganttData[ix].category;
                        const residueColor = ganttChart.get('colors').getIndex(ix)

                        selectSections.push({
                            residue_number: parseInt(residueID),
                            sideChain: true,
                            representation: 'spacefill',
                            representationColor: residueColor
                        })
                    }

                    viewerInstance.visual.select({
                        data: selectSections,
                    })
                });

            if (networkRootReference['active']) {
                lipidNodes = networkRootReference["series"].lipidNodes[lipidID]
                console.log('lipidNodes', lipidNodes)
                if (lipidNodes != undefined) {
                    am5.array.each(networkRootReference["series"].dataItems, function (dataItem) {
                        if (dataItem.dataContext.from != 0) {
                            dataItem.hide();
                        }
                        if (lipidNodes.includes(parseInt(dataItem.dataContext.from)) && lipidNodes.includes(parseInt(dataItem.dataContext.to))) {
                            dataItem.show();
                        }
                    });
                } else {
                    am5.array.each(networkRootReference["series"].dataItems, function (dataItem) {
                        if (dataItem.dataContext.from != 0) {
                            dataItem.hide();
                        }
                    });
                }
            }
        });

        // Should work for touch displays.
        table.on("rowTap", function (e, row) {
            obj = {
                "lipidID": row.getData()['lipidID']
            }
            fetch('/toplipids/' + JSON.stringify(obj))
                .then(response => response.json())
                .then(tableResponseData => {
                    ganttData = tableResponseData['ganttData'].map((lp, ix) => ({
                        ...lp
                    }))
                    ganttYAxis.data.setAll(tableResponseData['topLipids'].map(v => ({
                        category: v
                    })))
                    ganttSeries.data.setAll(ganttData);
                    sortCategoryAxis()
                });

        });

        // table.on("tableBuilt", function(e, row) {
        //     table.selectRow(0);
        // })

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
            wheelX: "none",
            wheelY: "none",
            layout: ganttRoot.verticalLayout
        }));
        ganttChart.zoomOutButton.set("forceHidden", true);

        ganttChart.children.unshift(am5.Label.new(ganttRoot, {
            text: "Lipid Contact Durations",
            x: am5.p50,
            centerX: am5.p50,
            centerY: am5.percent(25),
            position: "absolute",
            fontWeight: "500",
            fontStyle: "oblique",
        }));

        var legend = ganttChart.children.push(am5.Legend.new(ganttRoot, {
            centerX: am5.p50,
            x: am5.p50
        }))

        // var colors = ganttChart.get("colors");

        // This is kept as placeholder for setting custom color themes.
        // ganttData = responseData['ganttData'].map((lp, ix) => ({
        //     ...lp,
        //     columnSettings: {
        //         fill: colorSet.getIndex(ix * 3)
        //     }
        // }))
        ganttData = responseData['ganttData'].map((lp, ix) => ({
            ...lp
        }))

        var ganttYRenderer = am5xy.AxisRendererY.new(ganttRoot, {
            minGridDistance: 1,
            inversed: true,
        });

        // Create axes
        var ganttYAxis = ganttChart.yAxes.push(
            am5xy.CategoryAxis.new(ganttRoot, {
                categoryField: "category",
                maxDeviation: 0,
                renderer: ganttYRenderer,
            })
        );
        ganttYAxis.children.moveValue(am5.Label.new(ganttRoot, {
            text: "Residue IDs",
            rotation: -90,
            position: "absolute",
            // x: am5.p50,
            y: am5.percent(70),
            centerY: am5.percent(70),
            fontWeight: "500",
        }), 0);

        ganttYAxis.data.setAll(responseData['topLipids'].map(v => ({
            category: v,
        })))

        var ganttXAxis = ganttChart.xAxes.push(am5xy.ValueAxis.new(ganttRoot, {
            min: 0,
            max: frameNumber,
            // strictMinMax: true,
            renderer: am5xy.AxisRendererX.new(ganttRoot, {}),
            tooltip: am5.Tooltip.new(ganttRoot, {
                themeTags: ["axis"],
                animationDuration: 200
            })
        }));
        ganttXAxis.children.moveValue(am5.Label.new(ganttRoot, {
            text: "trajectory frames",
            position: "absolute",
            x: am5.p50,
            y: am5.p50,
            centerX: am5.p50,
            fontWeight: "500",
            // fontStyle: "oblique",
        }), ganttXAxis.children.length - 1);

        var ganttSeries = ganttChart.series.push(am5xy.ColumnSeries.new(ganttRoot, {
            xAxis: ganttXAxis,
            yAxis: ganttYAxis,
            openValueXField: "startFrame",
            valueXField: "endFrame",
            categoryYField: "category",
            sequencedInterpolation: true,
            tooltip: am5.Tooltip.new(ganttRoot, {
                pointerOrientation: "horizontal",
                labelText: "{category}"
            })

        }));

        ganttSeries.columns.template.setAll({
            // templateField: "columnSettings",
            strokeOpacity: 0,
            interactive: true,
            fillOpacity: 0.8,
            tooltipText: "{category}",
            cursorOverStyle: "pointer",
            // Rounded corners for bars
            cornerRadiusTR: 5,
            cornerRadiusBR: 5,
            cornerRadiusTL: 5,
            cornerRadiusBL: 5,

        });

        // Make each column to be of a different color
        ganttSeries.columns.template.adapters.add("fill", function (fill, target) {
            return ganttChart.get("colors").getIndex(ganttSeries.columns.indexOf(target));
        });

        ganttSeries.columns.template.adapters.add("stroke", function (stroke, target) {
            return ganttChart.get("colors").getIndex(ganttSeries.columns.indexOf(target));
        });

        ganttSeries.data.setAll(ganttData);
        sortCategoryAxis();

        csor = ganttChart.set("cursor", am5xy.XYCursor.new(ganttRoot, {
            behavior: "none",
            xAxis: ganttXAxis,
            yAxis: ganttYAxis
        }));

        // csor.lineY.set("visible", true);
        csor.lineX.set("opacity", 0.5);
        csor.lineY.set("opacity", 0.5);

        // Get series item by category
        function getSeriesItem(category) {
            for (var i = 0; i < ganttSeries.dataItems.length; i++) {
                var dataItem = ganttSeries.dataItems[i];
                if (dataItem.get("categoryY") == category) {
                    return dataItem;
                }
            }
        }

        // Axis sorting
        function sortCategoryAxis() {
            // Sort by value
            ganttSeries.dataItems.sort(function (x, y) {
                // return x.get("valueX") - y.get("valueX"); // descending
                return y.get("valueY") - x.get("valueX"); // ascending
            })

            // Go through each axis item
            am5.array.each(ganttYAxis.dataItems, function (dataItem) {
                // get corresponding series item
                var seriesDataItem = getSeriesItem(dataItem.get("category"));

                if (seriesDataItem) {
                    // get index of series data item
                    var index = ganttSeries.dataItems.indexOf(seriesDataItem);
                    // calculate delta position
                    var deltaPosition = (index - dataItem.get("index", 0)) / ganttSeries.dataItems.length;
                    // set index to be the same as series data item index
                    dataItem.set("index", index);
                    // set deltaPosition instanlty
                    dataItem.set("deltaPosition", -deltaPosition);
                    // animate delta position to 0
                    dataItem.animate({
                        key: "deltaPosition",
                        to: 0,
                        duration: 500,
                        easing: am5.ease.out(am5.ease.cubic)
                    })
                }
            });

            // Sort axis items by index.
            // This changes the order instantly, but as deltaPosition is set,
            // they keep in the same places and then animate to true positions.
            ganttYAxis.dataItems.sort(function (x, y) {
                return x.get("index") - y.get("index");
            });
        }

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

                    // Update title
                    var text = `ResidueID: ${ctx.category} and LipidID: ${ctx.lipid_id} Interactions`
                    am5.registry.entitiesById["besiTest"].set("text", text)

                    heatmapSeries.data.setAll(heatmapResponseData['heatmapData']);
                    hmYAxis.data.setAll(heatmapResponseData['residueAtomsData']);
                    hmXAxis.data.setAll(heatmapResponseData['lipidAtomsData']);
                });

            var cols = ganttChart.get('colors');
            var residueColor = am5.color("#E74C3C");
            am5.array.each(ganttSeries.dataItems, function (dataItem, ix) {
                if (dataItem.dataContext.category == residueID) {
                    residueColor = ganttChart.get('colors').getIndex(ix)
                }
            });

            var selectSections = [{
                residue_number: parseInt(ctx.category),
                // color:{r:255,g:0,b:255},
                representation: 'spacefill',
                representationColor: residueColor,
            }]
            viewerInstance.visual.select({
                data: selectSections,
            })
        });

        ganttSeries.columns.template.events.on("pointerover", function (e) {
            residueID = e.target.dataItem.dataContext.category;
            ctx = e.target.dataItem.dataContext;

            var selectSections = [{
                residue_number: parseInt(ctx.category),
                color: {
                    r: 255,
                    g: 0,
                    b: 255
                },
                representation: 'spacefill'
            }]
            viewerInstance.visual.highlight({
                data: selectSections,
            })
        })

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

        heatmapChart.children.unshift(am5.Label.new(heatmapRoot, {
            text: `ResidueID: ${startResidueID} and LipidID: ${startLipidID} Interactions`,
            id: "besiTest",
            x: am5.p50,
            centerX: am5.percent(40),
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




        // ##########################
        function networkApp(lipid = "CHOL") {

            var root = am5.Root.new("chartdiv");

            root.setThemes([
                am5themes_Animated.new(root)
            ]);

            var linkDefaultOpacity = 0.5,
                lihnkHoveredOpacity = 1;

            // Create series
            var series = root.container.children.push(am5flow.ChordNonRibbon.new(root, {
                sourceIdField: "from",
                targetIdField: "to",
                valueField: "value",
                padAngle: 0,
                startAngle: 90,
                draggable: true,
                lipid: lipid,

            }));

            series.nodes.labels.template.setAll({
                textType: "radial",
                fontSize: "0.6em",
                radius: 15,
            });

            series.nodes.bullets.push(function(_root, _series, dataItem) {
                var bulletCircle = am5.Circle.new(root, {
                  radius: 2.5,
                  fill: am5.color("#8E8A8A"),
                  fillOpacity: 1,
                });

                bulletCircle.adapters.add("fill", function(fill, target) {
                  var dataItem = target.dataItem;
                  if (dataItem) {
                    var sum = dataItem.get("sum", 0) // we can also use sumIncoming or sumOutgoing
                    var min = Infinity;
                    var max = -Infinity;
                    am5.array.each(series.nodes.dataItems, function(dataItem) {
                      var value = dataItem.get("sum");
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
                  sprite: bulletCircle
                });
              });

            series.children.moveValue(series.bulletsContainer, 0);

            hoverColor = am5.color(0xff0000);
            am5.array.each(subSeries.dataItems, function (dataItem, ix) {
                if (dataItem.dataContext.category == lipid) {
                    hoverColor = subSeries.get("colors").getIndex(ix)
                }
            })

            // Add heat legend
            var networkHeatLegend = series.children.push(am5.HeatLegend.new(root, {
                orientation: "horizontal",
                startColor: am5.color(0x673AB7),
                endColor: am5.color(0xF44336),
                width: 160,
                x: am5.percent(70),
                y: am5.percent(95),
                opacity: 1,
            }));

            series.nodes.rectangles.template.events.on("pointerover", function (ev) {
                var di = ev.target.dataItem;
                console.log('ev', di)
                if (di) {
                    networkHeatLegend.showValue(di._settings.sum);
                }

                var incomingLinks = ev.target.dataItem._settings.incomingLinks
                var outgoingLinks = ev.target.dataItem._settings.outgoingLinks
                if (incomingLinks != undefined) {
                    incomingLinks.forEach(link => {
                        link._settings.link._settings.stroke = hoverColor
                        link._settings.link._display.alpha = lihnkHoveredOpacity
                    })
                }
                if (outgoingLinks != undefined) {
                    outgoingLinks.forEach(link => {
                        link._settings.link._settings.stroke = hoverColor
                        link._settings.link._display.alpha = lihnkHoveredOpacity
                    })
                }
            })

            series.nodes.rectangles.template.events.on("pointerout", function (ev) {
                var incomingLinks = ev.target._dataItem._settings.incomingLinks
                var outgoingLinks = ev.target.dataItem._settings.outgoingLinks
                if (incomingLinks != undefined) {
                    incomingLinks.forEach(link => {
                        link._settings.link._settings.stroke = am5.color("#8E8A8A")
                        link._settings.link._display.alpha = linkDefaultOpacity
                    })
                }
                if (outgoingLinks != undefined) {
                    outgoingLinks.forEach(link => {
                        link._settings.link._settings.stroke = am5.color("#8E8A8A")
                        link._settings.link._display.alpha = linkDefaultOpacity
                    })
                }
            })

            obj = {
                "lipid": lipid,
                // "residueID": ctx.category
            }
            fetch('/network/' + JSON.stringify(obj))
                .then(response => response.json())
                .then(responseData => {

                    data = responseData['chordElements']
                    posRes = responseData['positionResidues']
                    series.data.setAll(data);
                    series.lipidNodes = responseData['lipidNodes'];

                    series.links.template.setAll({
                        strokeWidth: 0.2,
                        opacity: linkDefaultOpacity,
                        stroke: am5.color("#8E8A8A"),
                        strokeStyle: "none",
                    })

                    series.nodes.rectangles.template.setAll({
                        fillOpacity: 0,
                        fill: hoverColor,
                        tooltipText: "Residue [bold]{name}[/]\nShared Contacts: {sum}"
                    });

                    series.events.on("datavalidated", function () {

                        var bulletSums = [];
                        for (let ix = 0; ix < series.nodes.dataItems.length; ix++) {
                            var el = series.nodes.dataItems[ix];
                            sum = el.get('sum');
                            if (sum != 0) {
                                bulletSums.push(sum);
                            }
                        }

                        networkHeatLegend.set("startValue", Math.min(...bulletSums));
                        networkHeatLegend.set("endValue", Math.max(...bulletSums));

                        for (let ix = 0; ix < root.container.allChildren().length; ix++) {
                            const el = root.container.allChildren()[ix];

                            for (let jx = 0; jx < posRes.length; jx++) {
                                const pos = posRes[jx];
                                el.nodes.dataItems[pos].bullets[0]._settings.sprite._display.visible = false
                                el.nodes.labels._values[pos]._display.visible = false
                                el.nodes.rectangles.template._entities[pos]._display.visible = false
                            }
                            el.links._values.forEach((nodeLinks, ix) => {
                                if (nodeLinks.dataItem.dataContext.from == 0) {
                                    nodeLinks._settings.strokeWidth = 0;
                                    nodeLinks._settings.strokeOpacity = 0;
                                    nodeLinks._display.visible = false
                                } else {
                                    nodeLinks._settings.strokeWidth = nodeLinks.dataItem.dataContext.valueWidth * 2;
                                }
                            })
                        }
                    });
                });

            // Make stuff animate on load
            series.appear(1000, 100);


            var button2 = series.children.push(am5.Button.new(root, {
                dx: 0,
                dy: 0,
                label: am5.Label.new(root, {
                    text: "Circular Application",
                    fontSize: 10,
                })
            }));

            button2.events.on("click", function (ev) {
                root.dispose();
                fetch('/data/' + JSON.stringify({
                        lipid: lipid,
                        protein: "GIRK",
                    }))
                    .then(response => response.json())
                    .then(cirData => {

                        rootReferenceObjects = radarApp()

                        updateData = cirData['data']
                        rootReferenceObjects["series"].data.setAll(updateData);
                        rootReferenceObjects["categoryAxis"].data.setAll(updateData);
                        // rootReferenceObjects["createRange"]("POPE", updateData, 0);

                        am5.array.each(subSeries.dataItems, function (dataItem, ix) {
                            if (dataItem._settings.slice._settings.active) {
                                col = subSeries.get("colors").getIndex(ix)
                                rootReferenceObjects["createRange"](dataItem.dataContext.category, updateData, 0);
                                rootReferenceObjects["axisRange"].get("axisFill").set("fill", col)

                            }
                        })

                    });

                rootReferenceObjects["series"].appear(100);
            });

            return {
                "series": series,
                "root": root,
                "active": true
            }
        }

        ///////////////////////////////////////////
        ////////////// Hide Logos /////////////////
        ///////////////////////////////////////////
        am5.array.each(am5.registry.rootElements, function (rootElement) {
            rootElement.events.on("framestarted", function () {
                rootChildren = rootElement.tooltipContainer.allChildren()
                for (let ix = 0; ix < rootChildren.length; ix++) {
                    el = rootChildren[ix];
                    if (el._settings.tooltipText == "Created using amCharts 5") {
                        el.set('visible', false)
                    }
                }
            });
        });

    });