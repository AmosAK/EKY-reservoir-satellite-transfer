// =====================================================================
// EKY reservoirs: satellite-in situ matchups (Sentinel-2 + Landsat 8/9)
// Amos Akwemoh, CGIES/KSU
// For each station visit, extracts surface reflectance from every clear
// image within +/-3 days, as a 3x3-pixel median over water pixels only.
// =====================================================================

// 1. EDIT THIS: paste your asset path (Assets tab > click the table > copy ID)
var visits = ee.FeatureCollection('projects/YOUR_PROJECT/assets/eky_station_visits');

var WINDOW = 3;   // days either side of the sampling date
var MIN_PX = 5;   // minimum clear water pixels out of 9

print('Station visits:', visits.size());
Map.centerObject(visits, 8);
Map.addLayer(visits, {color: 'red'}, 'Stations');

// ---------- Sentinel-2 L2A (surface reflectance, 10-20 m) ----------
var S2_BANDS = ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'];

function prepS2(img) {
  var scl = img.select('SCL');
  // drop no-data, saturated, cloud shadow, medium/high cloud, cirrus
  var clear = scl.neq(0).and(scl.neq(1)).and(scl.neq(3))
                .and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  var water = img.normalizedDifference(['B3', 'B11']).gt(0);   // MNDWI
  return img.select(S2_BANDS).divide(10000)
            .updateMask(clear).updateMask(water)
            .set('sensor', img.get('SPACECRAFT_NAME'))
            .set('cloud_pct', img.get('CLOUDY_PIXEL_PERCENTAGE'))
            .copyProperties(img, ['system:time_start']);
}

// ---------- Landsat 8/9 Collection 2 Level-2 (30 m) ----------
var L_BANDS = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'];

function prepL(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)        // dilated cloud
               .and(qa.bitwiseAnd(1 << 3).eq(0))  // cloud
               .and(qa.bitwiseAnd(1 << 4).eq(0)); // cloud shadow
  var sr = img.select(L_BANDS).multiply(0.0000275).add(-0.2);
  var water = sr.normalizedDifference(['SR_B3', 'SR_B6']).gt(0);
  return sr.updateMask(clear).updateMask(water)
           .set('sensor', img.get('SPACECRAFT_ID'))
           .set('cloud_pct', img.get('CLOUD_COVER'))
           .copyProperties(img, ['system:time_start']);
}

// ---------- generic matchup extractor ----------
function matchups(collection, prepFn, scale, countBand) {
  return visits.map(function (v) {
    var d = ee.Date(v.get('date'));
    var box = v.geometry().buffer(scale * 1.5);   // ~3x3 pixels
    var imgs = collection.filterBounds(v.geometry())
                         .filterDate(d.advance(-WINDOW, 'day'), d.advance(WINDOW + 1, 'day'))
                         .map(prepFn);
    return imgs.map(function (img) {
      var med = img.reduceRegion({reducer: ee.Reducer.median(), geometry: box,
                                  scale: scale, maxPixels: 1e6});
      var n = img.select(countBand).reduceRegion({reducer: ee.Reducer.count(),
                                  geometry: box, scale: scale}).get(countBand);
      var t = ee.Date(img.get('system:time_start'));
      return ee.Feature(null, med)
               .set(v.toDictionary())
               .set({img_date: t.format('YYYY-MM-dd'),
                     dt_days: t.difference(d, 'day'),
                     n_px: n,
                     sensor: img.get('sensor'),
                     cloud_pct: img.get('cloud_pct')});
    });
  }).flatten().filter(ee.Filter.gte('n_px', MIN_PX));
}

var s2 = matchups(ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED'), prepS2, 10, 'B4');
var l89 = matchups(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                     .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')), prepL, 30, 'SR_B4');

print('Sentinel-2 matchups:', s2.size());
print('Landsat 8/9 matchups:', l89.size());
print('S2 example rows:', s2.limit(3));

// ---------- exports (Tasks tab > Run) ----------
Export.table.toDrive({collection: s2,  description: 'EKY_S2_matchups',
                      folder: 'EKY_reservoirs', fileFormat: 'CSV'});
Export.table.toDrive({collection: l89, description: 'EKY_L89_matchups',
                      folder: 'EKY_reservoirs', fileFormat: 'CSV'});
