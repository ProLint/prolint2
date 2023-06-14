import { sortCategoryAxis } from './ganttApp.js';
///////////////////////////////////////////
////////////// Lipid Table ////////////////
///////////////////////////////////////////
export function tableApp(responseData, ganttReturnValue, heatmap, networkRootReference) {

    var [ganttChart, ganttSeries, ganttXAxis, ganttYAxis] = ganttReturnValue;
    var [heatmapRoot, heatmapChart, heatmapSeries, hmYAxis, hmXAxis, heatLegend] = heatmap;

    var table = new Tabulator("#lipid-table", {
        data: responseData['tableData'],
        height: "290px",
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
                    formatter: function(cell) {
                        return parseFloat(cell.getValue()).toFixed(2);
                    },
                
                },
            ]
        }],
        headerVisible: true,
    });

    table.on("rowClick", function (e, row) {

        var classList = row.getElement().classList;

        var lipidID = row.getData()['lipidID']
        obj = {
            "lipidID": lipidID
        }
        if (classList.contains('tabulator-selected')) {
            fetch('/toplipids/' + JSON.stringify(obj))
                .then(response => response.json())
                .then(tableResponseData => {
                    var ganttData = tableResponseData['ganttData'].map((lp, ix) => ({
                        ...lp,
                    }))
                    ganttYAxis.data.setAll(tableResponseData['topLipids'].map(v => ({
                        category: v
                    })))
                    ganttSeries.data.setAll(ganttData);
                    sortCategoryAxis(ganttReturnValue)

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
                    obj = {
                        "lipidID": lipidID,
                        "residueID": ganttData[0].category
                    }
            
                    fetch('/distance/' + JSON.stringify(obj))
                    .then(response => response.json())
                    .then(heatmapResponseData => {
        
                        // Update title
                        var text = `Interactions between ResidueID: ${ganttData[0].category} and LipidID: ${lipidID}`
                        am5.registry.entitiesById["besiTest"].set("text", text)
        
                        heatmapSeries.data.setAll(heatmapResponseData['heatmapData']);
                        hmYAxis.data.setAll(heatmapResponseData['residueAtomsData']);
                        hmXAxis.data.setAll(heatmapResponseData['lipidAtomsData']);
                    });
        

                    viewerInstance.visual.select({
                        data: selectSections,
                    })
                });

        } else {
            viewerInstance.visual.select({
                data: [],
            })

        }

        if (document.getElementById("button_x").classList.contains("active")) {
        // if (!rootReferenceObjects['active']) {
            if (classList.contains('tabulator-selected')) {
                var lipidNodes = networkRootReference["series"].lipidNodes[lipidID]
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
            } else {
                am5.array.each(networkRootReference["series"].dataItems, function (dataItem) {
                    if (dataItem.dataContext.from != 0) {
                        dataItem.show();
                    }
                });
                // }
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
                var ganttData = tableResponseData['ganttData'].map((lp, ix) => ({
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

    return table;
}