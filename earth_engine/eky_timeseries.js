// =====================================================================
// EKY reservoirs: full Sentinel-2 time series, 2017-2026
// Amos Akwemoh, CGIES/KSU
// Output 1: station-level series (3x3 px, same footprint as the matchups,
//           so the calibrated models apply directly)
// Output 2: lake-wide series (median over the whole open-water area)
// Nothing to edit: just Run, then start both tasks in the Tasks tab.
// =====================================================================

var START = '2017-03-28', END = '2026-10-01';

// ---------- stations: [id, lake, lon, lat] ----------
var ST = [
  ['21KY_WQX-DOW04046004', 'Buckhorn', -83.461, 37.33023],
  ['21KY_WQX-DOW04046005', 'Buckhorn', -83.43671, 37.2744],
  ['21KY_WQX-DOW04046006', 'Buckhorn', -83.37058, 37.27204],
  ['21KY_WQX-CLN056', 'Cannon Creek', -83.70222, 36.68083],
  ['21KY_WQX-CLN322', 'Cannon Creek', -83.69266, 36.68478],
  ['21KY_WQX-CLN323', 'Cannon Creek', -83.7107, 36.67953],
  ['NARS_WQX-NLA_KY-10003', 'Cannon Creek', -83.70503, 36.68168],
  ['21KY_WQX-DOW05036025', 'Cave Run', -83.43877, 38.04373],
  ['21KY_WQX-DOW05036026', 'Cave Run', -83.51316, 38.09825],
  ['21KY_WQX-DOW05036029', 'Cave Run', -83.49677, 37.99698],
  ['21KY_WQX-DOW05036030', 'Cave Run', -83.4648, 38.05928],
  ['21KY_WQX-CLN055', 'Chenoa', -83.84944, 36.67583],
  ['21KY_WQX-CLN317', 'Chenoa', -83.84971, 36.67685],
  ['21KY_WQX-CLN057', 'Cranks Creek', -83.23758, 36.73907],
  ['21KY_WQX-CLN324', 'Cranks Creek', -83.22828, 36.74369],
  ['21KY_WQX-CLN325', 'Cranks Creek', -83.21686, 36.7476],
  ['21KY_WQX-CLN135', 'Dewey', -82.7277, 37.73563],
  ['21KY_WQX-CLN136', 'Dewey', -82.7216, 37.69842],
  ['21KY_WQX-CLN137', 'Dewey', -82.68197, 37.69831],
  ['21KY_WQX-CLN138', 'Fishtrap', -82.41504, 37.43408],
  ['21KY_WQX-CLN139', 'Fishtrap', -82.3687, 37.42486],
  ['21KY_WQX-CLN140', 'Fishtrap', -82.31185, 37.40982],
  ['21KY_WQX-CLN212', 'Grayson', -82.98513, 38.2522],
  ['21KY_WQX-CLN213', 'Grayson', -83.00806, 38.2339],
  ['21KY_WQX-CLN214', 'Grayson', -83.02999, 38.19633],
  ['21KY_WQX-DOW02043009', 'Martins Fork', -83.26104, 36.73735],
  ['21KY_WQX-CLN189', 'Paintsville', -82.872, 37.84101],
  ['21KY_WQX-CLN190', 'Paintsville', -82.926, 37.88299],
  ['21KY_WQX-CLN191', 'Paintsville', -82.988, 37.896],
  ['21KY_WQX-CLN192', 'Paintsville', -82.98839, 37.89822],
  ['NALMS-F509958', 'Paintsville', -82.889, 37.85],
  ['21KY_WQX-CLN062', 'Smoky Valley', -83.1183, 38.36629],
  ['21KY_WQX-CLN331', 'Smoky Valley', -83.12946, 38.3663],
  ['21KY_WQX-CLN193', 'Yatesville', -82.701, 38.124],
  ['21KY_WQX-CLN194', 'Yatesville', -82.743, 38.08901],
  ['21KY_WQX-CLN195', 'Yatesville', -82.78431, 38.05319],
  ['NALMS-F2742357', 'Yatesville', -82.772, 38.06],
];

var stations = ee.FeatureCollection(ST.map(function (s) {
  return ee.Feature(ee.Geometry.Point([s[2], s[3]]), {station: s[0], lake: s[1]});
}));

// ---------- lake polygons from JRC surface water (>=80% occurrence) ----------
// keep only water bodies that touch a station, then pull 20 m in from shore
var jrc = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').gte(80).selfMask();
var lakeNames = ['Buckhorn','Cannon Creek','Cave Run','Chenoa','Cranks Creek','Dewey',
                 'Fishtrap','Grayson','Martins Fork','Paintsville','Smoky Valley','Yatesville'];

var lakes = ee.FeatureCollection(lakeNames.map(function (name) {
  var pts = stations.filter(ee.Filter.eq('lake', name));
  var aoi = pts.geometry().buffer(9000).bounds();
  var vec = jrc.reduceToVectors({geometry: aoi, scale: 30, maxPixels: 1e10,
                                 geometryType: 'polygon', eightConnected: true});
  var poly = vec.filterBounds(pts.geometry().buffer(100)).geometry().buffer(-20, 10);
  var nTot = ee.Image.constant(1).reduceRegion({reducer: ee.Reducer.count(),
              geometry: poly, scale: 20, maxPixels: 1e10}).get('constant');
  return ee.Feature(poly, {lake: name, n_total: nTot,
                           area_km2: poly.area(10).divide(1e6)});
}));

Map.centerObject(stations, 8);
Map.addLayer(lakes, {color: '00bcd4'}, 'Lake polygons');
Map.addLayer(stations, {color: 'red'}, 'Stations');
print('Lakes (check area_km2 looks sensible):', lakes);

// ---------- Sentinel-2 L2A prep (same masks as the matchup script) ----------
var BANDS = ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'];
function prep(img) {
  var scl = img.select('SCL');
  var clear = scl.neq(0).and(scl.neq(1)).and(scl.neq(3))
                .and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  var water = img.normalizedDifference(['B3', 'B11']).gt(0);
  var swirOK = img.select('B11').divide(10000).lt(0.03);   // same QC as matchups
  return img.select(BANDS).divide(10000)
            .updateMask(clear).updateMask(water).updateMask(swirOK)
            .set({img_date: img.date().format('YYYY-MM-dd'),
                  sensor: img.get('SPACECRAFT_NAME'),
                  tile: img.get('MGRS_TILE'),
                  cloud_pct: img.get('CLOUDY_PIXEL_PERCENTAGE')})
            .copyProperties(img, ['system:time_start']);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
           .filterBounds(lakes.geometry())
           .filterDate(START, END)
           .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80))
           .map(prep);
print('Sentinel-2 scenes:', s2.size());

// ---------- Output 1: station-level series ----------
var stBuf = stations.map(function (f) { return f.buffer(15); });
var medCount = ee.Reducer.median().combine(ee.Reducer.count(), null, true);

var stationTS = s2.map(function (img) {
  return img.reduceRegions({collection: stBuf, reducer: medCount, scale: 10})
            .map(function (f) {
              return f.setGeometry(null).set({img_date: img.get('img_date'),
                sensor: img.get('sensor'), tile: img.get('tile'),
                cloud_pct: img.get('cloud_pct')});
            });
}).flatten().filter(ee.Filter.gte('B4_count', 5));

// ---------- Output 2: lake-wide series ----------
var lakeTS = lakes.map(function (lk) {
  return s2.filterBounds(lk.geometry()).map(function (img) {
    var st = img.reduceRegion({reducer: medCount, geometry: lk.geometry(),
                               scale: 20, maxPixels: 1e10, tileScale: 4});
    var n = ee.Number(st.get('B4_count'));
    return ee.Feature(null, st).set({lake: lk.get('lake'),
      img_date: img.get('img_date'), sensor: img.get('sensor'),
      tile: img.get('tile'), cloud_pct: img.get('cloud_pct'),
      area_km2: lk.get('area_km2'),
      frac_valid: n.divide(ee.Number(lk.get('n_total')))});
  });
}).flatten().filter(ee.Filter.gte('frac_valid', 0.3));

// ---------- exports (Tasks tab > Run both) ----------
Export.table.toDrive({collection: stationTS, description: 'EKY_S2_station_timeseries',
                      folder: 'EKY_reservoirs', fileFormat: 'CSV'});
Export.table.toDrive({collection: lakeTS, description: 'EKY_S2_lake_timeseries',
                      folder: 'EKY_reservoirs', fileFormat: 'CSV'});
Export.table.toDrive({collection: lakes, description: 'EKY_lake_polygons',
                      folder: 'EKY_reservoirs', fileFormat: 'GeoJSON'});
