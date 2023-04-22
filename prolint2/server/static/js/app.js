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
import { getTimeData  } from "./timeseries.js";
// NOTES (and TODO:)
// 1.
// The network app is always redrawn (rather than updated)
// This is mainly due to the pointerover events which retrieve the color
// value from the subseries pie chart. A solution might be update these
// events (or remove and add new events everytime the chart is updated).
// 2. Data fetching is not optimized to the scope of the different objects. The JSON.stringify
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
        var table = tableApp(responseData, ganttReturnValue, heatmap, networkRootReference)

        rootReferenceObjects["series"].data.setAll(contactData);
        rootReferenceObjects["categoryAxis"].data.setAll(contactData);
        rootReferenceObjects["createRange"](lipids[0], contactData, 0);

        rootReferenceObjects["series"].appear(100);
        rootReferenceObjects["chart"].appear(100);

        var [pieRoot, subSeries] = pieApp(table, ganttReturnValue, heatmap, timeSeries, networkRootReference, responseData, rootReferenceObjects);
   

        const selectElement = document.getElementById("metric_button");
        const selectedElementValue = selectElement.value;

        function selectSlice(slice) {
            const selectedSlice = slice;
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
        }
    
        function myFunction() {
          
            // Get the current active slice in the subSeries
            console.log('subSeries', subSeries)
            console.log('selectElement', document.getElementById("metric_button").value)
            // get lipid name from subseries
            subSeries.slices.values.forEach(function (slice) {

                if (slice.get('active')) {
                    console.log('slice', slice.dataItem.dataContext.category)


                    var lipid = slice.dataItem.dataContext.category;
                    obj.protein = "Protein";
                    obj.lipid = lipid;
                    obj.metric = document.getElementById("metric_button").value; 
                
                    fetch('/data/' + JSON.stringify(obj))
                        .then(response => response.json())
                        .then(responseData => {
                        var updateData = responseData['data'];
                        rootReferenceObjects["series"].data.setAll(updateData);
                        rootReferenceObjects["categoryAxis"].data.setAll(updateData);
                        rootReferenceObjects["createRange"](lipid, updateData, 0);
                
                        am5.array.each(subSeries.dataItems, function (dataItem, ix) {
                            if (dataItem.dataContext.category == lipid) {
                            var col = subSeries.get("colors").getIndex(ix);
                            rootReferenceObjects["axisRange"].get("axisFill").set("fill", col);
                            }
                        });
                
                        var [xAxis, series] = timeSeries;
                        var timeData = getTimeData(updateData);
                        xAxis.data.setAll(timeData);
                        series.data.setAll(timeData);
                        });
                    




                }
                // console.log('-> ', selectSlice(slice))
                // console.log('slice', slice)
                // console.log('slice', slice.dataItem.dataContext.category, slice.get('active'))
            })
            // var currentSlice = subSeries.slices.values.find(function (slice) {

            //     // console.log('-> ', selectSlice(slice))
            //     console.log('slice', slice.dataItem.dataContext.category)
            // //   return slice.isActive();
            // return true; 
            // });
          
            // if (currentSlice) {
            //   var lipid = currentSlice.dataItem.dataContext.category;
            //   obj.protein = "Protein";
            //   obj.lipid = lipid;
          
            //   fetch('/data/' + JSON.stringify(obj))
            //     .then(response => response.json())
            //     .then(responseData => {
            //       var updateData = responseData['data'];
            //       rootReferenceObjects["series"].data.setAll(updateData);
            //       rootReferenceObjects["categoryAxis"].data.setAll(updateData);
            //       rootReferenceObjects["createRange"](lipid, updateData, 0);
          
            //       am5.array.each(subSeries.dataItems, function (dataItem, ix) {
            //         if (dataItem.dataContext.category == lipid) {
            //           var col = subSeries.get("colors").getIndex(ix);
            //           rootReferenceObjects["axisRange"].get("axisFill").set("fill", col);
            //         }
            //       });
          
            //       var [xAxis, series] = timeSeries;
            //       var timeData = getTimeData(updateData);
            //       xAxis.data.setAll(timeData);
            //       series.data.setAll(timeData);
            //     });
            // }
          }

          selectElement.addEventListener("change", myFunction);

    

        ///////////////////////////////////////////
        ////////////// Hide Logos /////////////////
        ///////////////////////////////////////////
        const logos_to_keep_for_ids = ["chartdiv_x", "chartdiv_y", "chartdiv_z"]
        am5.array.each(am5.registry.rootElements, function (rootElement) {
            if (logos_to_keep_for_ids.includes(rootElement.dom.id)) {
                return
            }
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