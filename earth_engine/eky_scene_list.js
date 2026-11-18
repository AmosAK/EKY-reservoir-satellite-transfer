// =====================================================================
// List the exact satellite scenes behind the EKY test matchups
// (needed to download Level-1 data for ACOLITE).
// Output: one row per matchup with the scene product ID.
// =====================================================================
// [lon, lat, sensor, station, image date]
var ROWS = [
  [-83.23758, 36.73907, 'S2', '21KY_WQX-CLN057', '2023-07-10'],
  [-83.23758, 36.73907, 'S2', '21KY_WQX-CLN057', '2023-09-08'],
  [-83.1183, 38.36629, 'S2', '21KY_WQX-CLN062', '2023-07-23'],
  [-82.7216, 37.69842, 'S2', '21KY_WQX-CLN136', '2017-04-25'],
  [-82.68197, 37.69831, 'S2', '21KY_WQX-CLN137', '2017-04-25'],
  [-82.68197, 37.69831, 'S2', '21KY_WQX-CLN137', '2017-10-02'],
  [-82.872, 37.84101, 'S2', '21KY_WQX-CLN189', '2017-10-02'],
  [-82.98839, 37.89822, 'S2', '21KY_WQX-CLN192', '2017-05-02'],
  [-82.98839, 37.89822, 'S2', '21KY_WQX-CLN192', '2017-10-02'],
  [-82.701, 38.124, 'S2', '21KY_WQX-CLN193', '2017-10-17'],
  [-82.743, 38.08901, 'S2', '21KY_WQX-CLN194', '2017-10-17'],
  [-82.98513, 38.2522, 'S2', '21KY_WQX-CLN212', '2017-10-14'],
  [-83.00806, 38.2339, 'S2', '21KY_WQX-CLN213', '2017-10-17'],
  [-83.02999, 38.19633, 'S2', '21KY_WQX-CLN214', '2017-10-17'],
  [-83.22828, 36.74369, 'S2', '21KY_WQX-CLN324', '2023-07-10'],
  [-83.22828, 36.74369, 'S2', '21KY_WQX-CLN324', '2023-09-06'],
  [-83.21686, 36.7476, 'S2', '21KY_WQX-CLN325', '2023-07-10'],
  [-83.21686, 36.7476, 'S2', '21KY_WQX-CLN325', '2023-09-06'],
  [-83.12946, 38.3663, 'S2', '21KY_WQX-CLN331', '2023-07-23'],
  [-82.772, 38.06, 'S2', 'NALMS-F2742357', '2019-07-24'],
  [-82.772, 38.06, 'S2', 'NALMS-F2742357', '2019-08-05'],
  [-82.772, 38.06, 'S2', 'NALMS-F2742357', '2019-08-20'],
  [-82.772, 38.06, 'S2', 'NALMS-F2742357', '2019-09-09'],
  [-82.772, 38.06, 'S2', 'NALMS-F2742357', '2019-10-24'],
  [-82.772, 38.06, 'S2', 'NALMS-F2742357', '2020-06-10'],
  [-82.889, 37.85, 'S2', 'NALMS-F509958', '2019-06-21'],
  [-82.889, 37.85, 'S2', 'NALMS-F509958', '2019-07-09'],
  [-82.889, 37.85, 'S2', 'NALMS-F509958', '2019-07-24'],
  [-82.889, 37.85, 'S2', 'NALMS-F509958', '2019-08-05'],
  [-83.70222, 36.68083, 'L', '21KY_WQX-CLN056', '2023-07-18'],
  [-83.69266, 36.68478, 'L', '21KY_WQX-CLN322', '2023-07-18'],
  [-83.7107, 36.67953, 'L', '21KY_WQX-CLN323', '2023-07-18'],
  [-83.70222, 36.68083, 'L', '21KY_WQX-CLN056', '2023-09-04'],
  [-83.69266, 36.68478, 'L', '21KY_WQX-CLN322', '2023-09-04'],
  [-83.7107, 36.67953, 'L', '21KY_WQX-CLN323', '2023-09-04'],
  [-83.84944, 36.67583, 'L', '21KY_WQX-CLN055', '2023-09-04'],
  [-83.23758, 36.73907, 'L', '21KY_WQX-CLN057', '2023-07-10'],
  [-83.21686, 36.7476, 'L', '21KY_WQX-CLN325', '2023-07-11'],
  [-83.22828, 36.74369, 'L', '21KY_WQX-CLN324', '2023-09-04'],
  [-83.21686, 36.7476, 'L', '21KY_WQX-CLN325', '2023-09-04'],
  [-82.7277, 37.73563, 'L', '21KY_WQX-CLN135', '2017-07-17'],
  [-82.7216, 37.69842, 'L', '21KY_WQX-CLN136', '2017-07-17'],
  [-82.68197, 37.69831, 'L', '21KY_WQX-CLN137', '2017-07-17'],
  [-82.7277, 37.73563, 'L', '21KY_WQX-CLN135', '2017-10-05'],
  [-82.7216, 37.69842, 'L', '21KY_WQX-CLN136', '2017-10-05'],
  [-82.68197, 37.69831, 'L', '21KY_WQX-CLN137', '2017-10-05'],
  [-82.41504, 37.43408, 'L', '21KY_WQX-CLN138', '2017-07-17'],
  [-83.02999, 38.19633, 'L', '21KY_WQX-CLN214', '2017-07-17'],
  [-82.872, 37.84101, 'L', '21KY_WQX-CLN189', '2017-10-05'],
  [-83.12946, 38.3663, 'L', '21KY_WQX-CLN331', '2023-07-26'],
  [-82.701, 38.124, 'L', '21KY_WQX-CLN193', '2017-07-10'],
  [-82.743, 38.08901, 'L', '21KY_WQX-CLN194', '2017-10-14'],
  [-82.772, 38.06, 'L', 'NALMS-F2742357', '2019-08-08'],
  [-82.772, 38.06, 'L', 'NALMS-F2742357', '2019-08-24'],
  [-82.772, 38.06, 'L', 'NALMS-F2742357', '2019-09-09'],
  [-82.772, 38.06, 'L', 'NALMS-F2742357', '2019-10-11'],
  [-82.772, 38.06, 'L', 'NALMS-F2742357', '2020-06-07']
];
var pts = ee.FeatureCollection(ee.List(ROWS).map(function (r) {
  r = ee.List(r);
  return ee.Feature(ee.Geometry.Point(ee.List([r.get(0), r.get(1)])),
                    {sensor: r.get(2), station: r.get(3), img_date: r.get(4)});
}));

function sceneFor(f) {
  var d = ee.Date(f.get('img_date'));
  var isS2 = ee.String(f.get('sensor')).equals('S2');
  var s2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED').filterBounds(f.geometry()).filterDate(d, d.advance(1, 'day'));
  var ls = ee.ImageCollection('LANDSAT/LC08/C02/T1').merge(ee.ImageCollection('LANDSAT/LC09/C02/T1'))
             .filterBounds(f.geometry()).filterDate(d, d.advance(1, 'day'));
  var ids = ee.Algorithms.If(isS2, s2.aggregate_array('PRODUCT_ID'), ls.aggregate_array('LANDSAT_PRODUCT_ID'));
  return f.set('product_ids', ee.List(ids).join(';'));
}
var out = pts.map(sceneFor);
print('Scenes:', out.limit(10));
Export.table.toDrive({collection: out, description: 'EKY_test_scene_list', folder: 'EKY_reservoirs_HLS',
                      fileFormat: 'CSV', selectors: ['sensor', 'station', 'img_date', 'product_ids']});
