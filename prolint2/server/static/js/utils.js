export function sortCategoryAxis2(ganttReturnValue) {

    var [ganttChart, ganttSeries, ganttXAxis, ganttYAxis] = ganttReturnValue;
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
