// 在 GEE Code Editor 中运行。
// 导出 2024 年 1--5 月、10--12 月逐小时 GSMaP 站点长表。

var YEAR = 2024;
var MONTHS = [1, 2, 3, 4, 5, 10, 11, 12];
var START = ee.Date.fromYMD(YEAR, 1, 1);
var END = START.advance(1, 'year'); // 右开区间

var AOI = ee.Geometry.Rectangle(
  [108.672, 29.232, 116.029, 33.200], null, false);

var stationsRaw = ee.FeatureCollection(
  'projects/booming-list-489917-p4/assets/Hubei-STATIONS');

// 兼容站点资产中可能存在的 station_id/site/id 字段命名差异。
function pickKeyContains(feature, keyword) {
  var keys = ee.List(feature.propertyNames());
  var matched = keys.map(function(key) {
    key = ee.String(key);
    return ee.Algorithms.If(
      key.replace('_', '', 'g').toLowerCase().index(keyword).gte(0), key, null);
  }).removeAll([null]);
  return ee.Algorithms.If(
    matched.size().gt(0), feature.get(ee.String(matched.get(0))), null);
}

function normalizeStation(feature) {
  var hiddenSite = pickKeyContains(feature, 'site');
  var stationId = ee.Algorithms.If(
    feature.get('station_id'), feature.get('station_id'),
    ee.Algorithms.If(
      feature.get('site'), feature.get('site'),
      ee.Algorithms.If(
        hiddenSite, hiddenSite,
        ee.Algorithms.If(
          feature.get('id'), feature.get('id'), feature.get('system:index')))));
  var lon = ee.Algorithms.If(
    feature.get('lon'), feature.get('lon'), pickKeyContains(feature, 'lon'));
  var lat = ee.Algorithms.If(
    feature.get('lat'), feature.get('lat'), pickKeyContains(feature, 'lat'));

  return ee.Feature(
    ee.Geometry.Point([ee.Number(lon), ee.Number(lat)]), {
      station_id: ee.String(stationId),
      lon: ee.Number(lon),
      lat: ee.Number(lat)
    });
}

// 不按 AOI 执行 filterBounds，避免浮点误差漏掉边界站点。
// 不对 station_id 去重，以保持项目原有月文件的记录结构。
var stations = stationsRaw
  .map(normalizeStation)
  .filter(ee.Filter.notNull(['station_id', 'lon', 'lat']));

// GSMaP V8 逐小时非雨量站订正版，hourlyPrecipRate 单位为 mm/hr。
var gsmap = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v8/operational')
  .filterDate(START, END)
  .select('hourlyPrecipRate');

function hourlyLongTable(monthStart, monthEnd) {
  var hourCount = monthEnd.difference(monthStart, 'hour');
  var hours = ee.List.sequence(0, hourCount.subtract(1));

  return ee.FeatureCollection(hours.map(function(offset) {
    var hourStart = monthStart.advance(ee.Number(offset), 'hour');
    var hourEnd = hourStart.advance(1, 'hour');
    // GSMaP 每小时一景；mean 保持该小时 mm/hr 数值，并兼容集合筛选结果。
    var image = gsmap.filterDate(hourStart, hourEnd).mean();
    var sampled = image.reduceRegions({
      collection: stations,
      reducer: ee.Reducer.mean(),
      scale: 11132,
      tileScale: 4
    });

    return sampled.map(function(feature) {
      return ee.Feature(feature.geometry(), {
        gsmap_mm_h: feature.get('mean'),
        station_id: feature.get('station_id'),
        time: hourStart.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
      });
    }).filter(ee.Filter.notNull(['station_id', 'gsmap_mm_h']));
  })).flatten();
}

function exportMonth(year, month) {
  var monthStart = ee.Date.fromYMD(year, month, 1);
  var monthEnd = monthStart.advance(1, 'month');
  var table = hourlyLongTable(monthStart, monthEnd);
  var suffix = String(year) + (month < 10 ? '0' + month : String(month));

  print('已创建导出任务：gsmap_hubei_hourly_long_' + suffix);
  Export.table.toDrive({
    collection: table,
    description: 'gsmap_hubei_hourly_long_' + suffix,
    folder: 'GEE_GSMAP_HUBEI',
    fileNamePrefix: 'gsmap_hubei_hourly_long_' + suffix,
    fileFormat: 'CSV',
    selectors: ['system:index', 'gsmap_mm_h', 'station_id', 'time', '.geo']
  });
}

MONTHS.forEach(function(month) {
  exportMonth(YEAR, month);
});

Map.centerObject(AOI, 7);
Map.addLayer(stations, {color: 'yellow'}, 'Hubei stations');
