var svg = document.childNodes[0];
var concat = String.prototype.concat ;
var sqrt2 = Math.sqrt(2) ;
var rec = function(alpha, x, y, theta, l, r, opp) {
  if (r < 0.5) { return; }
  var newx = x + Math.cos(theta)*l
  var newy = y - Math.sin(theta)*l
  var lr = 0.6

  var path = document.createElementNS("http://www.w3.org/2000/svg", "path");
  var str = "M ";
  str = str.concat(x, " ", y, " L ", newx," ", newy);
  path.setAttribute("d", str);
  path.setAttribute("stroke", "black");
  path.setAttribute("stroke-width", r);
  path.setAttribute("fill", "none");
  svg.appendChild(path);

  var maxasym = 4 ;
  var runif = Math.random() ;
  var asym = 1/(runif*(1 - 1/maxasym) + 1/maxasym)
  var minr = (1/(1 + asym**alpha))**(1/alpha)
  var majr = asym*minr
  var minl = (1/(1 + asym**2))**(1/2)
  var majl = asym*minl
  var mintheta = Math.atan(asym)
  var majtheta = Math.atan(1/asym)
  var rot = Math.cos(Math.random() * Math.PI)
  var wind = 0.9*Math.PI/2
  var windamt = 0.3
  var majtheta = wind*windamt + (1 - windamt)*(theta + rot*majtheta)
  var mintheta = wind*windamt + (1 - windamt)*(theta - rot*mintheta)
  rec(alpha, newx, newy, majtheta, l * majl, r*majr, !opp);
  rec(alpha, newx, newy, mintheta, l * minl, r*minr, !opp);
}
rec(1, 150, 400, Math.PI/2, 50, 5, true);
rec(2, 450, 400, Math.PI/2, 50, 5, true);
rec(3, 750, 400, Math.PI/2, 50, 5, true);
