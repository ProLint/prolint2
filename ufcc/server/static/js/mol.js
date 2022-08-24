viewerInstance.events.loadComplete.subscribe(() => {
    viewerInstance.canvas.toggleControls(false)
    viewerInstance.canvas.setBgColor({r:255, g:255, b:255})
    viewerInstance.plugin.canvas3d.camera.state.mode = "orthographic"
    // console.log('vI', viewerInstance)
    // viewerInstance.visual.select({
    //   data: [{
    //     entity_id: '*'
    //   }],
    //   nonSelectedColor: {r:255,g:0,b:0}
    // })
  });
// viewerInstance.canvas.setBgColor({r:255, g:255, b:255})
