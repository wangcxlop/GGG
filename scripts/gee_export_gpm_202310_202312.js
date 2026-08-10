// 在 GEE Code Editor 中运行。
// 作用：按项目现有的 GPM 长表格式，导出 2023 年 10、11、12 月的逐小时站点数据。
// 导出的 CSV 文件名：gpm_hubei_hourly_long_202310.csv（以及 202311、202312）。

var START = ee.Date('2023-10-01T00:00:00');
var END = ee.Date('2024-01-01T00:00:00'); // 右开区间

// 沿用项目中已确认的站点经纬度范围。
var AOI = ee.Geometry.Rectangle(
  [108.672, 29.232, 116.029, 33.200], null, false);

// 现有 Hubei GPM 脚本使用的站点资产。
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

// Hubei-STATIONS 已经是目标站点集合，因此不再按 AOI 过滤站点。
// 若按与站点极值完全相同的边界 filterBounds，浮点坐标误差会漏掉边界站点。
// 不对 station_id 做 distinct，保留现有 GPM 导出中站点资产的重复记录行为；
// 项目内 prepare_satellite_inputs.jl 会对相同 station_id/time 的相同值取均值。
var stations = stationsRaw
  .map(normalizeStation)
  .filter(ee.Filter.notNull(['station_id', 'lon', 'lat']));

// IMERG V07 为半小时产品；precipitation 单位为 mm/hr。
var gpm = ee.ImageCollection('NASA/GPM_L3/IMERG_V07')
  .filterDate(START, END)
  .select('precipitation');

function hourlyLongTable(monthStart, monthEnd) {
  var hourCount = monthEnd.difference(monthStart, 'hour');
  var hours = ee.List.sequence(0, hourCount.subtract(1));

  return ee.FeatureCollection(hours.map(function(offset) {
    var hourStart = monthStart.advance(ee.Number(offset), 'hour');
    var hourEnd = hourStart.advance(1, 'hour');
    // 两个半小时影像取均值，与现有 hourly_long 文件的 mm/hr 语义一致。
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

exportMonth(2023, 10, 2023, 11);
exportMonth(2023, 11, 2023, 12);
exportMonth(2023, 12, 2024, 1);

Map.centerObject(AOI, 7);
Map.addLayer(stations, {color: 'yellow'}, 'Hubei stations');
