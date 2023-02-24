// import { sortCategoryAxis } from './utils.js';
///////////////////////////////////////////
/////////////// GanttApp //////////////////
///////////////////////////////////////////
export function ganttApp(responseData, heatmap) {

    var frameNumber = responseData["frameNumber"];
    var [heatmapRoot, heatmapChart, heatmapSeries, hmYAxis, hmXAxis, heatLegend] = heatmap;

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

    var ganttData = responseData['ganttData'].map((lp, ix) => ({
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

    var ganttReturnValue = [ganttChart, ganttSeries, ganttXAxis, ganttYAxis];
    ganttSeries.data.setAll(ganttData);
    sortCategoryAxis(ganttReturnValue);

    var csor = ganttChart.set("cursor", am5xy.XYCursor.new(ganttRoot, {
        behavior: "none",
        xAxis: ganttXAxis,
        yAxis: ganttYAxis
    }));

    // csor.lineY.set("visible", true);
    csor.lineX.set("opacity", 0.5);
    csor.lineY.set("opacity", 0.5);


    ganttSeries.appear();
    ganttChart.appear(1000, 100);

    ganttSeries.columns.template.events.on("click", function (e, d) {
        var residueID = e.target.dataItem.dataContext.category;

        var ctx = e.target.dataItem.dataContext;

        obj = {
            "lipidID": ctx.lipid_id,
            "residueID": ctx.category
        }
        fetch('/distance/' + JSON.stringify(obj))
            .then(response => response.json())
            .then(heatmapResponseData => {

                // Update title
                var text = `Interactions between ResidueID: ${ctx.category} and LipidID: ${ctx.lipid_id}`
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
        var residueID = e.target.dataItem.dataContext.category;
        var ctx = e.target.dataItem.dataContext;

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

    return [ganttChart, ganttSeries, ganttXAxis, ganttYAxis]

}

export function sortCategoryAxis(ganttReturnValue) {
    
    var [ganttChart, ganttSeries, ganttXAxis, ganttYAxis] = ganttReturnValue;
    
    // Get series item by category
    function getSeriesItem(category) {
        for (var i = 0; i < ganttSeries.dataItems.length; i++) {
            var dataItem = ganttSeries.dataItems[i];
            if (dataItem.get("categoryY") == category) {
                return dataItem;
            }
        }
    }

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
