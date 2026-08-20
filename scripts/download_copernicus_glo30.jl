#!/usr/bin/env julia

using Dates
using Downloads
using SHA

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_OUTPUT = joinpath(ROOT, "data", "raw", "copernicus_glo30")

const S3_HOST = "eodata.dataspace.copernicus.eu"
const S3_ENDPOINT = "https://$S3_HOST"
const S3_BUCKET = "eodata"
const S3_REGION = "default"
const S3_SERVICE = "s3"
const S3_PREFIX = "auxdata/CopDEM_COG/copernicus-dem-30m"
const EMPTY_SHA256 = bytes2hex(sha256(UInt8[]))

const LATITUDES = 31:33
const LONGITUDES = 109:111

struct Config
    output::String
    overwrite::Bool
    dry_run::Bool
end

function usage(io::IO=stdout)
    println(io, """
    下载研究区所需的 9 个 Copernicus DEM GLO-30 官方 COG 瓦片。

    用法：
      julia --project=. scripts/download_copernicus_glo30.jl [选项]

    选项：
      --output DIR   下载目录（默认：data/raw/copernicus_glo30）
      --overwrite    覆盖已经存在的有效 TIFF 文件
      --dry-run      只列出瓦片和官方 S3 路径，不进行下载
      -h, --help     显示帮助

    下载前需要设置 CDSE 生成的 S3 凭据：
      CDSE_S3_ACCESS_KEY
      CDSE_S3_SECRET_KEY
    """)
end

function parse_args(args::Vector{String})
    output = DEFAULT_OUTPUT
    overwrite = false
    dry_run = false
    i = 1

    while i <= length(args)
        arg = args[i]
        if arg == "--output"
            i == length(args) && error("--output 后必须提供目录")
            i += 1
            output = abspath(args[i])
        elseif startswith(arg, "--output=")
            output = abspath(split(arg, "="; limit=2)[2])
        elseif arg == "--overwrite"
            overwrite = true
        elseif arg == "--dry-run"
            dry_run = true
        elseif arg == "-h" || arg == "--help"
            usage()
            return nothing
        else
            error("未知参数：$(arg)（使用 --help 查看帮助）")
        end
        i += 1
    end

    return Config(output, overwrite, dry_run)
end

tile_id(lat::Integer, lon::Integer) =
    "Copernicus_DSM_COG_10_N$(lpad(lat, 2, '0'))_00_E$(lpad(lon, 3, '0'))_00_DEM"

function tile_assets()
    return [
        let id = tile_id(lat, lon)
            (id=id, key="$S3_PREFIX/$id/$id.tif")
        end
        for lat in LATITUDES for lon in LONGITUDES
    ]
end

function aws_headers(access_key::String, secret_key::String, canonical_uri::String)
    timestamp = now(UTC)
    amz_date = Dates.format(timestamp, "yyyymmddTHHMMSSZ")
    date_stamp = Dates.format(timestamp, "yyyymmdd")
    signed_headers = "host;x-amz-content-sha256;x-amz-date"

    canonical_headers =
        "host:$S3_HOST\n" *
        "x-amz-content-sha256:$EMPTY_SHA256\n" *
        "x-amz-date:$amz_date\n"
    canonical_request = join(
        [
            "GET",
            canonical_uri,
            "",
            canonical_headers,
            signed_headers,
            EMPTY_SHA256,
        ],
        '\n',
    )

    scope = "$date_stamp/$S3_REGION/$S3_SERVICE/aws4_request"
    string_to_sign = join(
        [
            "AWS4-HMAC-SHA256",
            amz_date,
            scope,
            bytes2hex(sha256(canonical_request)),
        ],
        '\n',
    )

    date_key = hmac_sha256(Vector{UInt8}(codeunits("AWS4$secret_key")), date_stamp)
    region_key = hmac_sha256(date_key, S3_REGION)
    service_key = hmac_sha256(region_key, S3_SERVICE)
    signing_key = hmac_sha256(service_key, "aws4_request")
    signature = bytes2hex(hmac_sha256(signing_key, string_to_sign))

    authorization =
        "AWS4-HMAC-SHA256 Credential=$access_key/$scope, " *
        "SignedHeaders=$signed_headers, Signature=$signature"

    return [
        "Authorization" => authorization,
        "x-amz-content-sha256" => EMPTY_SHA256,
        "x-amz-date" => amz_date,
    ]
end

function is_tiff(path::AbstractString)
    isfile(path) || return false
    filesize(path) >= 4 || return false
    magic = open(path, "r") do io
        read(io, 4)
    end
    return magic in (
        UInt8[0x49, 0x49, 0x2a, 0x00], # little-endian TIFF
        UInt8[0x4d, 0x4d, 0x00, 0x2a], # big-endian TIFF
        UInt8[0x49, 0x49, 0x2b, 0x00], # little-endian BigTIFF
        UInt8[0x4d, 0x4d, 0x00, 0x2b], # big-endian BigTIFF
    )
end

function download_tile(
    asset,
    output_dir::AbstractString,
    access_key::String,
    secret_key::String;
    overwrite::Bool=false,
)
    destination = joinpath(output_dir, "$(asset.id).tif")
    partial = "$destination.part"

    if isfile(destination) && !overwrite
        is_tiff(destination) || error(
            "目标文件已存在但不是有效 TIFF：$(destination)。" *
            "请人工检查，确认后使用 --overwrite。",
        )
        println("  已存在，跳过：$(basename(destination))")
        return :skipped
    end

    isfile(partial) && rm(partial; force=true)
    canonical_uri = "/$S3_BUCKET/$(asset.key)"
    url = "$S3_ENDPOINT$canonical_uri"
    headers = aws_headers(access_key, secret_key, canonical_uri)

    response = Downloads.request(url; headers=headers, output=partial, throw=false)
    if !(200 <= response.status < 300)
        isfile(partial) && rm(partial; force=true)
        if response.status == 403
            error(
                "S3 返回 HTTP 403。请检查 CDSE S3 凭据，并确认账户已取得 " *
                "Copernicus DEM GLO-30 的访问许可。",
            )
        elseif response.status == 404
            error("S3 返回 HTTP 404，官方资产不存在：$(asset.key)")
        else
            error("下载失败，HTTP $(response.status)：$(asset.key)")
        end
    end

    is_tiff(partial) || error("下载结果不是有效 TIFF：$(partial)")
    mv(partial, destination; force=overwrite)
    mib = round(filesize(destination) / 1024^2; digits=1)
    println("  完成：$(basename(destination)) ($(mib) MiB)")
    return :downloaded
end

function credentials()
    access_key = String(strip(get(ENV, "CDSE_S3_ACCESS_KEY", "")))
    secret_key = String(strip(get(ENV, "CDSE_S3_SECRET_KEY", "")))
    isempty(access_key) && error("缺少环境变量 CDSE_S3_ACCESS_KEY")
    isempty(secret_key) && error("缺少环境变量 CDSE_S3_SECRET_KEY")
    return access_key, secret_key
end

function main(args::Vector{String}=ARGS)
    config = parse_args(args)
    config === nothing && return nothing

    assets = tile_assets()
    @assert length(assets) == 9

    println("研究区：109.4°E–111.6°E，31.2°N–33.4°N")
    println("瓦片范围：N31–N33、E109–E111，共 $(length(assets)) 个")
    println("保存目录：$(config.output)")

    if config.dry_run
        for (i, asset) in enumerate(assets)
            println("[$i/$(length(assets))] s3://$S3_BUCKET/$(asset.key)")
        end
        return assets
    end

    access_key, secret_key = credentials()
    mkpath(config.output)
    downloaded = 0
    skipped = 0

    for (i, asset) in enumerate(assets)
        println("[$i/$(length(assets))] $(asset.id)")
        flush(stdout)
        status = download_tile(
            asset,
            config.output,
            access_key,
            secret_key;
            overwrite=config.overwrite,
        )
        downloaded += status == :downloaded
        skipped += status == :skipped
    end

    println(
        "下载完成：新增 $(downloaded) 个，跳过 $(skipped) 个，总计 $(length(assets)) 个。",
    )
    println("原始 COG 保存在：$(config.output)")
    return assets
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
