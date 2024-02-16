var lengths = Array.prototype.slice.call(document.getElementsByTagName("path")).map(function(a){return a.getTotalLength()}).join("\n");
