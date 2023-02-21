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

import { networkApp } from "../js/network.js";
import { radarApp } from "../js/radialApp.js";
import { ganttApp } from "./ganttApp.js";
import { heatmapApp } from "./heatmapApp.js";
import { pieApp } from "./pie.js";
import { tableApp } from "./table.js";
import { timeSeriesApp } from "./timeseries.js";
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

// Fetch the data from the backend
var obj = {
    "lipid": "",
    "protein": ""
}

fetch('/data/' + JSON.stringify(obj))
    .then(response => response.json())
    .then(responseData => {

        console.log('responseData', responseData)

        var contactData = responseData['data'];
        var lipids = responseData['lipids'];

        var rootReferenceObjects = radarApp()

        var heatmap = heatmapApp(responseData);
        var timeSeries = timeSeriesApp(contactData); 
        var ganttReturnValue = ganttApp(responseData, heatmap);
        var networkRootReference = networkApp(subSeries, lipids[0])
        var table = tableApp(responseData, ganttReturnValue, networkRootReference)

        rootReferenceObjects["series"].data.setAll(contactData);
        rootReferenceObjects["categoryAxis"].data.setAll(contactData);
        rootReferenceObjects["createRange"](lipids[0], contactData, 0);

        rootReferenceObjects["series"].appear(100);
        rootReferenceObjects["chart"].appear(100);

        var [pieRoot, subSeries] = pieApp(table, ganttReturnValue, heatmap, timeSeries, networkRootReference, responseData, rootReferenceObjects);
   
        ///////////////////////////////////////////
        ////////////// Hide Logos /////////////////
        ///////////////////////////////////////////
        am5.array.each(am5.registry.rootElements, function (rootElement) {
            rootElement.events.on("framestarted", function () {
                var rootChildren = rootElement.tooltipContainer.allChildren()
                for (let ix = 0; ix < rootChildren.length; ix++) {
                    var el = rootChildren[ix];
                    if (el._settings.tooltipText == "Created using amCharts 5") {
                        el.set('visible', false)
                    }
                }
            });
        });

    });