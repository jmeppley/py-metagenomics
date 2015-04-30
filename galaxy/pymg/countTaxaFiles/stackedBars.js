var dataFile=getDataFileFromForm();
var transpose=getTranspose();
var dataFormat=getDataFileFormat();
var dataMaps = {}, yScales = {};
var layout = d3.select("[id=layoutSelect]").property("value");
var sortMethod = d3.select("[id=sortSelect]").property("value");

var fullWidth=1200;
var legendWidth=250;
var margin = {top: 20, right: 20, bottom: 30, left: 40},
    width = fullWidth - legendWidth - margin.left - margin.right - 1,
    height = 500 - margin.top - margin.bottom;

var xo = d3.scale.ordinal()
    .rangeRoundBands([0, width], .1);
var xl = d3.scale.linear()
        .range([0,width]);

var y = d3.scale.linear()
    .rangeRound([height, 0]);

var legendScale = d3.scale.ordinal();
var nameKey;

var color = d3.scale.category20();

// Set up a container to put SVG and legend into
var container = d3.select("body").append("div")
    .attr("id","container")
    .style("min-width",fullWidth + "px")
    .style("float","left")
    .style("overflow","hidden");

var svgWidth=width + margin.left + margin.right;
var svg = container.append("div")
    .style("float","left")
    .style("width",svgWidth+"px").append("svg")
    .attr("width",svgWidth) 
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

var legend = container.append("div")
    .attr("id","legend")
    .style("float","left")
    .style("width",legendWidth+"px");
    legend.append("div").style("height",margin.top+"px");

loadData(dataFormat, dataFile);

/* This is borrowed from d3.dsv.parse. I had to postopone appliaction 
    until after the data was transposed. */
function keydata(data) {
    var kdata=Array();
    var o;
    data.forEach(function(d,i) {
        if (o) {
            // This will only be tru AFTER the first/header line
            kdata.push(o(d));
        } else {
            /* returns something like "return {header1:d[0],header2:d[1]}" */
            o = new Function("d", "return {" + d.map(function(name, i) {
                        return JSON.stringify(name) + ": d[" + i + "]";
                              }).join(",") + "}");
        }

    });
    return kdata;
}

function getSortedLegendNames() {
    return d3.range(legendScale.range().length)
    .map(function(i) { 
        return legendScale(i); 
    });
}

var rowSum = function(row) {
    return d3.sum(d3.values(row));
};

function setLegendScaleDomain(sortMethod) {
    var names = legendScale.range().slice(0);
    var data = dataMaps;
  if (sortMethod != "none") {
      if (sortMethod.substring(sortMethod.length-3)=='pct') {
          valueFunction = new Function("data,name","return d3." + sortMethod.substring(0,sortMethod.length-3) + "(data.map(function(row) { return row[name]/rowSum(row); }));");
      } else {
          valueFunction = new Function("data,name","return d3." + sortMethod + "(data.map(function(row) { return row[name]; }));");
      }
      names.sort(function(a,b) { b=valueFunction(data,b); a=valueFunction(data,a); return b - a; });
  }
  /* return list suitable for ordinal scale domain
  * ie. 3rd element should be index that 3rd name should have in sorted list
  */
  //d3.select("body").append("p").text(JSON.stringify(names));
  legendScale.domain(legendScale.range().map(function(name) {
      return names.indexOf(name);
  }));
}

function calculateHeightsFull(row) {
    if (layout == 'fractional') {
        var yScale=d3.scale.linear()
          .domain([0,row.total])
          .rangeRound([height,0]);
    } else {
         var yScale=y;
    }

    // Initialize variable here to keep track of where stacked bars go
    var y0 = 0;
    // Generate array of positions in order given
    return d3.range(legendScale.range().length).map(function(index) {
        var name = legendScale(index);
        //links back to data label and parent row
        ys={name: name, p:row};
        // For each rectangle, set y0 to top of last, and height to val/tot
        ys['y0'] = yScale(y0);
        y0=y0*1.0+row[name]/1.0;
        ys['y1'] = yScale(y0);
        return ys;
    });
}

function setBarWidths() {
    if (layout=='fractional') {
      // Add all row totals
      var grandTotal = d3.sum(dataMaps.map(function(row) { return row.total; }));
      //d3.select("body").append("p").text(JSON.stringify(grandTotal));
      // set the x scale to map from 0 to grandTotal
      xl.domain([0,grandTotal*1.1]);
      var padding=0.1*grandTotal/dataMaps.length;
      var nextX=padding*0.5;
      dataMaps.forEach(function(d) {
        d.position=xl(nextX);
        d.width=xl(d.total);
        nextX = nextX + d.total + padding;
      });
    } else {
        // Use the range bands
        xo.domain(dataMaps.map(function(d) { return d[nameKey]; }));
        dataMaps.forEach(function(d) {
            d.position=xo(d[nameKey]);
            d.width=xo.rangeBand()
        });
        if (layout=='stacked') {
            y.domain([0,d3.max(dataMaps.map(function(row) { return row.total; }))]);
        }
   }
}

function createPlot(data) {

  nameKey=data[0][0];

  data = transpose ? d3.transpose(data) : data;
  var legendNames = data[0].slice(1);
  legendScale.range(legendNames);
  data = keydata(data);
  dataMaps = data;
  setLegendScaleDomain(sortMethod);
  color.domain(getSortedLegendNames());
  //d3.select("body").append("p").text(JSON.stringify(legendScale.range()));
  //d3.select("body").append("p").text(JSON.stringify(legendScale.domain()));

  // Create a 'total' for each row
  data.forEach(function(row) { row.total=d3.sum(legendNames.map(function(n) { return row[n]; })); });
  
  setBarWidths();

      // Create a container for each bar/stack
  var state = svg.selectAll(".state")
      .data(data)
    .enter().append("g")
      .attr("class", "state")
      ;

      // Create a rectangle for each element in a stacked bar
  state.selectAll("rect")
  .data(function(d) { return calculateHeightsFull(d); }, 
                      function(d) { return d.name; })
    .enter().append("rect")
    .attr("x", function(d) { return d.p.position; })
    .attr("width", function(d) { return d.p.width; })
         .attr("y", function(d) { return d.y1; })
         .attr("height", function(d) { return d.y0 - d.y1; })
      .style("fill", function(d) { return color(d.name); })
      .append("svg:title").text(function(d) { return d.name + ": " + d.p[d.name]; })
      ;

  state.append("text")
  .attr("class", "axis")
  .attr("transform", function(d) { return "translate(" + (d.width+d.position-5.0) + ","+(height-5)+") rotate(-90)"; })
      .text(function(d) { return d[nameKey]; });

  state.append("text")
   .attr("class", "data")
   .attr("text-anchor", "middle")
   .attr("transform", function(d) { return "translate(" + (d.position+d.width*0.5) + ")"; })
   .style("opacity","0")
   .text(function(d) { return JSON.stringify(d.total); })
   ;

legend.selectAll(".legend")
.data(getSortedLegendNames().reverse())
    .enter().append("div")
    .attr("class","legend")
    .style("color", function(name) { return color(name); })
    .text(function(name) { return name; })
    ;
}

function loadData(dataFormat, dataFile) {
    d3.select("[id=downloadLink]").attr("href",dataFile);
    if (dataFormat == "csv") {
        d3.csv.parse = d3.csv.parseRows;
        d3.csv(dataFile, function(error, data) {
            createPlot(data);
        });
    } else {
        d3.tsv.parse = d3.tsv.parseRows;
        d3.tsv(dataFile, function(error, data) {
            createPlot(data);
        });
    }
}

function clearData() {
    svg.selectAll(".state").remove();
    legend.selectAll(".legend").remove();
}

function changeData() {
    dataFile=getDataFileFromForm();
    clearData();
    loadData(dataFormat, dataFile);
}

function changeLayout() {
    layout=getLayoutFromForm();
    setBarWidths();

    // apply changes
    var state = svg.selectAll(".state");

    state.selectAll("rect")
     .data(function(d) { return calculateHeightsFull(d); },
                      function(d) { return d.name; })
     .transition()
        .duration(750)
        .attr("x", function(d) { return d.p.position; })
        .attr("width", function(d) { return d.p.width; })
     .transition()
         .duration(750)
         .attr("y", function(d) { return d.y1; })
         .attr("height", function(d) { return d.y0 - d.y1; });
    state
     .selectAll(".axis")
      .transition(750)
       .attr("transform", function(d) { return "translate(" + (d.width+d.position-5.0) + ","+(height-5)+") rotate(-90)"; })
        
     ;
    state
     .selectAll(".data")
      .transition(1750)
       .style("opacity",layout=="stacked" ? "1" : "0")
      .transition(1750)
       .attr("transform", function(d) { return "translate(" + (d.position+d.width*0.5) + "," + ((layout=="stacked") ? y(d.total)-5 : 0) + ")"; })
      
} 

function changeSort() {
    var sortMethod = getSortMethodFromForm();
    setLegendScaleDomain(sortMethod);

    // Re-order rectangles in bars
    svg.selectAll(".state")
     .selectAll("rect")
     .data(function(d) { return calculateHeightsFull(d); },
                      function(d) { return d.name; })
     .transition()
     .duration(1000)
     .attr("y", function(d) { return d.y1; })
     .attr("height", function(d) { return d.y0 - d.y1; })
     ;

     // re-sort legend (w/out animation for now)
    legend.selectAll(".legend")
    .data(getSortedLegendNames().reverse())
        .style("color", function(name) { return color(name); })
        .text(function(name) { return name; })
        ;

}


