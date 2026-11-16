// =====================================================================
// HLS matchups: Harmonized Landsat Sentinel-2 (HLS v2.0) for the same
// station visits used in the original analysis.
// HLSS30 (Sentinel-2) and HLSL30 (Landsat 8/9) share one atmospheric
// correction (LaSRC), a common 30 m grid, BRDF normalisation and
// bandpass adjustment, so differences between them are sensor-only.
// Same rules as the original scripts: +/-3 days, 3x3-pixel median,
// MNDWI > 0 water pixels, >= 5 valid pixels. dt_days is exported so the
// window can be narrowed later (same day, +/-1 day) without re-running.
// =====================================================================

// ---- 1. Asset paths (same tables as the original matchup scripts) ----
var PROJECT = 'projects/extreme-clone-482209-b7/assets/';
var visitsReg = ee.FeatureCollection(PROJECT + 'regional_station_visits_s2era');
var visitsNLA = ee.FeatureCollection(PROJECT + 'nla_2017_2022_visits');
var visitsEKY = ee.FeatureCollection(PROJECT + 'eky_station_visits');   // check this name in your Assets tab

var WINDOW = 3;   // days either side of sampling date
var MIN_PX = 5;   // minimum valid water pixels out of 9
var SCALE = 30;   // HLS grid (m)

// ---- 2. Reflectance scaling: CHECKPOINT ----
// Run the script once and look at the printed 'HLS sample values'.
// If red/green values look like 0.01-0.2, set REFL_SCALE = 1.
// If they look like 100-2000, set REFL_SCALE = 0.0001.
var REFL_SCALE = 1;

// ---- 3. Common band names for both sensors ----
// HLSS30 (Sentinel-2): B1 coastal, B2 blue, B3 green, B4 red, B5-B7 red edge,
//                      B8A narrow NIR, B11 SWIR1, B12 SWIR2
// HLSL30 (Landsat):    B1 coastal, B2 blue, B3 green, B4 red, B5 NIR,
//                      B6 SWIR1, B7 SWIR2
var S30_IN  = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12'];
var S30_OUT = ['coastal', 'blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'swir1', 'swir2'];
var L30_IN  = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7'];
var L30_OUT = ['coastal', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2'];

// Fmask bits: 1 cloud, 2 adjacent to cloud/shadow, 3 cloud shadow, 4 snow/ice
function clearMask(img) {
  var f = img.select('Fmask');
  return f.bitwiseAnd(1 << 1).eq(0)
          .and(f.bitwiseAnd(1 << 2).eq(0))
          .and(f.bitwiseAnd(1 << 3).eq(0))
          .and(f.bitwiseAnd(1 << 4).eq(0));
}

function prep(inBands, outBands, sensorName) {
  return function (img) {
    var sr = img.select(inBands, outBands).multiply(REFL_SCALE);
    var water = sr.normalizedDifference(['green', 'swir1']).gt(0);   // MNDWI
    return sr.updateMask(clearMask(img)).updateMask(water)
             .set('sensor', sensorName)
             .set('cloud_pct', img.get('CLOUD_COVERAGE'))
             .copyProperties(img, ['system:time_start']);
  };
}

var S30 = ee.ImageCollection('NASA/HLS/HLSS30/v002');
var L30 = ee.ImageCollection('NASA/HLS/HLSL30/v002');

// ---- 4. Matchup extractor (identical logic to the original scripts) ----
function matchups(visits, collection, prepFn) {
  return visits.map(function (v) {
    var d = ee.Date(v.get('date'));
    var box = v.geometry().buffer(SCALE * 1.5);   // ~3x3 pixels
    var imgs = collection.filterBounds(v.geometry())
                         .filterDate(d.advance(-WINDOW, 'day'), d.advance(WINDOW + 1, 'day'))
                         .map(prepFn);
    return imgs.map(function (img) {
      var med = img.reduceRegion({reducer: ee.Reducer.median(), geometry: box,
                                  scale: SCALE, maxPixels: 1e6});
      var n = img.select('red').reduceRegion({reducer: ee.Reducer.count(),
                                  geometry: box, scale: SCALE}).get('red');
      var t = ee.Date(img.get('system:time_start'));
      return ee.Feature(null, med)
               .set(v.toDictionary())
               .set({img_date: t.format('YYYY-MM-dd'),
                     img_time_utc: t.format('HH:mm'),
                     dt_days: t.difference(d, 'day'),
                     n_px: n,
                     sensor: img.get('sensor'),
                     cloud_pct: img.get('cloud_pct')});
    });
  }).flatten().filter(ee.Filter.gte('n_px', MIN_PX));
}

var prepS30 = prep(S30_IN, S30_OUT, 'HLS_S30');
var prepL30 = prep(L30_IN, L30_OUT, 'HLS_L30');

// ---- 5. Checkpoint print (small sample only) ----
var sample = matchups(visitsEKY.limit(20), S30, prepS30).limit(3);
print('HLS sample values (check red/green scale):', sample);

// ---- 6. Exports (Tasks tab > Run each) ----
var sets = [['REG', visitsReg], ['NLA', visitsNLA], ['EKY', visitsEKY]];
sets.forEach(function (s) {
  Export.table.toDrive({collection: matchups(s[1], S30, prepS30),
                        description: 'HLS_S30_' + s[0] + '_matchups',
                        folder: 'EKY_reservoirs_HLS', fileFormat: 'CSV'});
  Export.table.toDrive({collection: matchups(s[1], L30, prepL30),
                        description: 'HLS_L30_' + s[0] + '_matchups',
                        folder: 'EKY_reservoirs_HLS', fileFormat: 'CSV'});
});
