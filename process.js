const fs = require('fs');
const os = require('os');
const path = require('path');

const RGRAPH = /^Loading graph .*\/(.+?)\.mtx \.\.\./m;
const RORDER = /^order: (\d+) size: (\d+) \[directed\] \{\} \(transposeWithDegree\)$/m;
const ROMPTH = /^OMP_NUM_THREADS=(\d+)/m;
const RLOOPF = /^pagerank\w+LoopU\(\):/m;
const RTSTAT = /^\[thread (\d+)\] status {processed=(\d+), stolen=(\d+), slept=(\d+), blocked=(\S+) ms}/m;
const RRESLT = /^\[(.+?) ms; (.+?) iters\.\] \[(.+?) err\.\] (\w+)(?:\s+\{sleep_prob: (.+?), sleep_dur: (\d+) ms\})?/m;




// UTILITY
// -------

function arraySum(x) {
  var a = 0;
  for (var v of x)
    a += v;
  return a;
}




// *-FILE
// ------

function readFile(pth) {
  var d = fs.readFileSync(pth, 'utf8');
  return d.replace(/\r?\n/g, '\n');
}

function writeFile(pth, d) {
  d = d.replace(/\r?\n/g, os.EOL);
  fs.writeFileSync(pth, d);
}




// *-CSV
// -----

function writeCsv(pth, rows) {
  var cols = Object.keys(rows[0]);
  var a = cols.join()+'\n';
  for (var r of rows)
    a += [...Object.values(r)].map(v => `"${v}"`).join()+'\n';
  writeFile(pth, a);
}




// *-LOG
// -----

function blankRowCount(datum) {
  var a = 0;
  for (var i=datum.length-1; i>=0; --i) {
    if (datum[i].technique) break;
    ++a;
  }
  return a;
}

function correctedTime(datum, time, repeat) {
  var B = blankRowCount(datum), db = Math.ceil(B / repeat), maxs = [];
  for (var i=datum.length-B; i<datum.length; i+=db)
    maxs.push(Math.max(...datum.slice(i, i+db).map(x => x.blocked_time)));
  return time - arraySum(maxs)/maxs.length;
}

function readLogLine(ln, data, state) {
  if (RGRAPH.test(ln)) {
    var [, graph] = RGRAPH.exec(ln);
    if (!data.has(graph)) data.set(graph, []);
    state = {graph};
  }
  else if (RORDER.test(ln)) {
    var [, order, size] = RORDER.exec(ln);
    state.order = parseFloat(order);
    state.size  = parseFloat(size);
  }
  else if (ROMPTH.test(ln)) {
    var [, omp_num_threads] = ROMPTH.exec(ln);
    state.omp_num_threads   = parseFloat(omp_num_threads);
    state.repeat            = 0;
  }
  else if (RLOOPF.test(ln)) {
    state.repeat++;
    state.time       = 0;
    state.iterations = 0;
    state.error      = 0;
    state.technique  = '';
    state.sleep_probability = 0;
    state.sleep_duration    = 0;
  }
  else if (RTSTAT.test(ln)) {
    var [, thread, processed_count, stolen_count, slept_count, blocked_time] = RTSTAT.exec(ln);
    data.get(state.graph).push(Object.assign({}, state, {
      thread: parseFloat(thread),
      processed_count: parseFloat(processed_count),
      stolen_count:    parseFloat(stolen_count),
      slept_count:     parseFloat(slept_count),
      blocked_time:    parseFloat(blocked_time),
    }));
  }
  else if (RRESLT.test(ln)) {
    var [, time, iterations, error, technique, sleep_probability, sleep_duration] = RRESLT.exec(ln);
    var result = {
      time:       parseFloat(time),
      iterations: parseFloat(iterations),
      error:      parseFloat(error),
      technique,
      sleep_probability: parseFloat(sleep_probability || '0'),
      sleep_duration:    parseFloat(sleep_duration    || '0'),
    };
    if (!result.technique.includes('Barrierfree')) result.corrected_time = result.time;
    else result.corrected_time = correctedTime(data.get(state.graph), result.time, state.repeat);
    var datum = data.get(state.graph);
    for (var i=datum.length-1; i>=0; --i) {
      if (datum[i].technique) break;
      Object.assign(datum[i], result);
    }
    state.repeat = 0;
    state.thread = 0;
    state.processed_count = 0;
    state.stolen_count    = 0;
    state.slept_count     = 0;
    state.blocked_time    = 0;
    data.get(state.graph).push(Object.assign({}, state, result));
  }
  return state;
}

function readLog(pth) {
  var text = readFile(pth);
  var lines = text.split('\n');
  var data = new Map();
  var state = null;
  for (var ln of lines)
    state = readLogLine(ln, data, state);
  return data;
}




// PROCESS-*
// ---------

function processCsv(data) {
  var a = [];
  for (var rows of data.values())
    a.push(...rows);
  return a;
}




// MAIN
// ----

function main(cmd, log, out) {
  var data = readLog(log);
  if (path.extname(out)==='') cmd += '-dir';
  switch (cmd) {
    case 'csv':
      var rows = processCsv(data);
      writeCsv(out, rows);
      break;
    case 'csv-dir':
      for (var [graph, rows] of data)
        writeCsv(path.join(out, graph+'.csv'), rows);
      break;
    default:
      console.error(`error: "${cmd}"?`);
      break;
  }
}
main(...process.argv.slice(2));
