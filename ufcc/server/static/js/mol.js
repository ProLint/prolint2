viewerInstance.events.loadComplete.subscribe(() => {
    viewerInstance.canvas.toggleControls(false)
    viewerInstance.canvas.setBgColor({r:255, g:255, b:255})
    viewerInstance.plugin.canvas3d.camera.state.mode = "orthographic"
  });

// get toggle button
var controltoggle = document.querySelectorAll('[title="Toggle Controls Panel"]')[0];
// change background color
if (controltoggle.parentElement) {
  controltoggle.style.backgroundColor = "#f0ebeb";
}
// lower opacity of child element (svg drawing):
controltoggle.childNodes[0].style.opacity = 0.6;
// disable pointer events:
controltoggle.style.pointerEvents = "none";

var fullScreenToggle = document.querySelectorAll('[title="Toggle Expanded Viewport"]')[0];
fullScreenToggle.addEventListener('click', function(ev) {
  if (Array.from(fullScreenToggle.classList).includes("msp-btn-link-toggle-off")) {
    controltoggle.style.backgroundColor = "transparent";
    controltoggle.childNodes[0].style.opacity = 1;
    controltoggle.style.pointerEvents = "auto";
  } else {
    controltoggle.style.backgroundColor = "#f0ebeb";
    controltoggle.childNodes[0].style.opacity = 0.6;
    controltoggle.style.pointerEvents = "none";

  }
})
