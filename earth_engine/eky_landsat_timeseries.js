// =====================================================================
// EKY reservoirs: Landsat 8/9 station time series, 2017-2026
// Same stations, dates and footprint logic as the Sentinel-2 series,
// so the two sensors can be compared over the exact same years.
// Amos Akwemoh, CGIES/KSU
// Nothing to edit: Run, then start the task in the Tasks tab.
// Output: EKY_L89_station_timeseries.csv in Drive folder EKY_reservoirs
// =====================================================================

var START = '2017-03-28', END = '2026-10-01';   // identical to eky_timeseries.js

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
  ['NALMS-F2742357', 'Yatesville', -82.772, 38.06]
];
var stations = ee.FeatureCollection(ST.map(function (s) {
  return ee.Feature(ee.Geometry.Point([s[2], s[3]]), {station: s[0], lake: s[1]});
}));

// ---------- Landsat 8/9 Collection 2 L2 prep (same rules as the matchup script) ----------
var SR  = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'];
var OUT = ['coastal','blue','green','red','nir','swir1','swir2'];
function prep(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)        // dilated cloud
               .and(qa.bitwiseAnd(1 << 3).eq(0)) // cloud
               .and(qa.bitwiseAnd(1 << 4).eq(0)) // cloud shadow
               .and(qa.bitwiseAnd(1 << 5).eq(0));// snow
  var sr = img.select(SR).multiply(0.0000275).add(-0.2).rename(OUT);
  var water  = sr.normalizedDifference(['green', 'swir1']).gt(0);   // MNDWI > 0
  var swirOK = sr.select('swir1').lte(0.03);                        // same QC as training
  return sr.updateMask(clear).updateMask(water).updateMask(swirOK)
           .set({img_date: img.date().format('YYYY-MM-dd'),
                 sensor: img.get('SPACECRAFT_ID'),
                 path_row: ee.String(img.get('WRS_PATH')).cat('_').cat(ee.String(img.get('WRS_ROW'))),
                 cloud_pct: img.get('CLOUD_COVER')})
           .copyProperties(img, ['system:time_start']);
}
var l89 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
            .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
            .filterBounds(stations.geometry().buffer(5000))
            .filterDate(START, END)
            .filter(ee.Filter.lt('CLOUD_COVER', 80))
            .map(prep);
print('Landsat 8/9 scenes:', l89.size());

// ---------- station series: 3x3 pixels at 30 m (45 m buffer), >=5 valid water pixels ----------
var stBuf = stations.map(function (f) { return f.buffer(45); });
var medCount = ee.Reducer.median().combine(ee.Reducer.count(), null, true);
var ts = l89.map(function (img) {
  return img.reduceRegions({collection: stBuf, reducer: medCount, scale: 30, tileScale: 2})
            .map(function (f) {
              return f.setGeometry(null).set({img_date: img.get('img_date'),
                sensor: img.get('sensor'), path_row: img.get('path_row'),
                cloud_pct: img.get('cloud_pct')});
            });
}).flatten().filter(ee.Filter.gte('red_count', 5));

print('Example rows:', ts.limit(5));
Map.centerObject(stations, 8);
Map.addLayer(stations, {color: 'red'}, 'Stations');

Export.table.toDrive({
  collection: ts, description: 'EKY_L89_station_timeseries',
  folder: 'EKY_reservoirs', fileFormat: 'CSV',
  selectors: ['station','lake','img_date','sensor','path_row','cloud_pct',
              'coastal_median','blue_median','green_median','red_median',
              'nir_median','swir1_median','swir2_median','red_count']
});
