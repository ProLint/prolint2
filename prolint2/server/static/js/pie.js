import { networkApp } from "../js/network.js";
import { getTimeData  } from "./timeseries.js";
export function pieApp(table, ganttReturnValue, heatmap, timeSeries, networkRootReference, responseData, rootReferenceObjects) {


    var [ganttChart, ganttSeries, ganttXAxis, ganttYAxis] = ganttReturnValue; 
    var [heatmapRoot, heatmapChart, heatmapSeries, hmYAxis, hmXAxis, heatLegend] = heatmap;

    var lipids = responseData['lipids'];
    var proteins = responseData['proteins'];

    var systemHasOneProtein = false;
    if (proteins.length == 1) {
        systemHasOneProtein = true
    }

    
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

                var updateData = responseData['data']
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
        // if (document.getElementById("button_y").classList.contains("active")) {
        if (lipid != rootReferenceObjects["axisRange"].get("label").get('text')) {
            sameLipid = false;
            // TODO:
            // get correct protein
            // Update Circular App Data
            obj.protein = "Protein"
            obj.lipid = lipid
            obj.metric = document.getElementById("metric_button").value;
            fetch('/metric/' + JSON.stringify(obj))
                .then(response => response.json())
                .then(responseData => {

                    var updateData = responseData['data'];

                    rootReferenceObjects["series"].data.setAll(updateData);
                    rootReferenceObjects["categoryAxis"].data.setAll(updateData);
                    rootReferenceObjects["createRange"](lipid, updateData, 0);

                    am5.array.each(subSeries.dataItems, function (dataItem, ix) {
                        if (dataItem.dataContext.category == lipid) {
                            var col = subSeries.get("colors").getIndex(ix)
                            rootReferenceObjects["axisRange"].get("axisFill").set("fill", col)
                        }
                    })

                    var [xAxis, series] = timeSeries
                    var timeData = getTimeData(updateData);
                    xAxis.data.setAll(timeData);
                    series.data.setAll(timeData);
                
                });
        }

        sameLipid = false;
        networkRootReference["root"].dispose();
        networkRootReference = networkApp(lipid)
        networkRootReference["series"].appear(1000)

        // } 

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

                            var ganttData = tableResponseData['ganttData'].map((lp, ix) => ({
                                ...lp
                            }))
                            ganttYAxis.data.setAll(tableResponseData['topLipids'].map(v => ({
                                category: v
                            })))
                            ganttSeries.data.setAll(ganttData);
                            // sortCategoryAxis()

                            // Update Heatmap App Data
                            obj = {
                                "lipidID": lipid_id,
                                "residueID": ganttData[0]['category']
                            }
            
                            fetch('/distance/' + JSON.stringify(obj))
                                .then(response => response.json())
                                .then(heatmapResponseData => {

                                    var text = `ResidueID: ${ganttData[0]['category']} and LipidID: ${lipid_id} Interactions`
                                    am5.registry.entitiesById["besiTest"].set("text", text)
        
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

    return [pieRoot, subSeries]; 

}