/* javascript funtions for turning on/off details boxes
 * Mario Stanke
 */
function onoff (linkDivID) {
  if (!document.getElementById)
     return;
  var linkItem = document.getElementById(linkDivID);
  var showItem = document.getElementById(linkItem.title);

  if (!showItem || !linkItem)
     return;
  if (showItem.style.display=='none')
     turnOn(showItem, linkItem);
  else
     turnOff(showItem, linkItem);
}

function turnOff (showItem, linkItem) {
     showItem.style.display='none';
     linkItem.innerHTML = "[+]";
}
function turnOn (showItem, linkItem) {
     showItem.style.display='block';
     linkItem.innerHTML = "[-]";
}
function allOn(){
     var details = document.getElementsByClassName('dcross');
     for (var i = 0; (d = details[i]) != null; i++) {
         turnOn(document.getElementById(d.title), d);
     }
}
function allOff(){
     var details = document.getElementsByClassName('dcross');
     for (var i = 0; (d = details[i]) != null; i++) {
         turnOff(document.getElementById(d.title), d);
     }
}
