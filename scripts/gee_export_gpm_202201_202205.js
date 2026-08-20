// 在 GEE Code Editor 中运行。
// 导出 2022 年 1--5 月逐小时 GPM 站点长表。

var START = ee.Date('2022-01-01T00:00:00');
var END = ee.Date('2022-06-01T00:00:00'); // 右开区间
var AOI = ee.Geometry.Rectangle(
  [109.4, 31.2, 111.6, 33.4], null, false);

var stationsRaw = ee.FeatureCollection(
  'projects/booming-list-489917-p4/assets/Hubei-STATIONS');

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

// 使用经纬度数值闭区间筛选，避免几何边界浮点误差。
var stations = stationsRaw
  .map(normalizeStation)
  .filter(ee.Filter.notNull(['station_id', 'lon', 'lat']))
  .filter(ee.Filter.gte('lon', 109.4))
  .filter(ee.Filter.lte('lon', 111.6))
  .filter(ee.Filter.gte('lat', 31.2))
  .filter(ee.Filter.lte('lat', 33.4));

var gpm = ee.ImageCollection('NASA/GPM_L3/IMERG_V07')
  .filterDate(START, END)
  .select('precipitation');

function hourlyLongTable(monthStart, monthEnd) {
  var hourCount = monthEnd.difference(monthStart, 'hour');
  var hours = ee.List.sequence(0, hourCount.subtract(1));

  return ee.FeatureCollection(hours.map(function(offset) {
    var hourStart = monthStart.advance(ee.Number(offset), 'hour');
    var hourEnd = hourStart.advance(1, 'hour');
    var image = gpm.filterDate(hourStart, hourEnd).mean();
    var sampled = image.reduceRegions({
      collection: stations,
      reducer: ee.Reducer.mean(),
      scale: 11132,
      tileScale: 4
    });

    return sampled.map(function(feature) {
      return ee.Feature(feature.geometry(), {
        gpm_mm_h: feature.get('mean'),
        station_id: feature.get('station_id'),
        time: hourStart.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
      });
    }).filter(ee.Filter.notNull(['station_id', 'gpm_mm_h']));
  })).flatten();
}

function exportMonth(year, month, nextYear, nextMonth) {
  var monthStart = ee.Date.fromYMD(year, month, 1);
  var monthEnd = ee.Date.fromYMD(nextYear, nextMonth, 1);
  var table = hourlyLongTable(monthStart, monthEnd);
  var suffix = String(year) + (month < 10 ? '0' + month : String(month));

  print('已创建导出任务：gpm_hubei_hourly_long_' + suffix);
  Export.table.toDrive({
    collection: table,
    description: 'gpm_hubei_hourly_long_' + suffix,
    folder: 'GEE_GPM_HUBEI',
    fileNamePrefix: 'gpm_hubei_hourly_long_' + suffix,
    fileFormat: 'CSV',
    selectors: ['system:index', 'gpm_mm_h', 'station_id', 'time', '.geo']
  });
}

exportMonth(2022, 1, 2022, 2);
exportMonth(2022, 2, 2022, 3);
exportMonth(2022, 3, 2022, 4);
exportMonth(2022, 4, 2022, 5);
exportMonth(2022, 5, 2022, 6);

Map.centerObject(AOI, 7);
Map.addLayer(stations, {color: 'yellow'}, 'Hubei stations');
