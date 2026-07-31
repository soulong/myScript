# ============================================================
# WSL Ubuntu 镜像定位 + 空间释放与压缩 + 备份 (PowerShell 7)
# 需以【管理员身份】运行
# 每步确认: y=执行  n=取消退出  s=跳过此步进入下一步
# ============================================================
$ErrorActionPreference = 'Stop'

# --- 确认辅助函数：返回 'y' / 'n' / 's' ---
function Confirm-Step {
    param(
        [string]$Message,
        [string]$Default = 'n'
    )
    while ($true) {
        $answer = Read-Host "$Message [y=执行 / n=取消 / s=跳过]"
        if ([string]::IsNullOrWhiteSpace($answer)) { $answer = $Default }
        switch ($answer.ToLower()) {
            'y' { return 'y' }
            'n' { return 'n' }
            's' { return 's' }
            default { Write-Host "请输入 y / n / s。" }
        }
    }
}

# --- 1. 找到 Ubuntu 发行版名称 ---
$raw    = wsl --list --quiet
$distro = $raw |
    ForEach-Object { ($_ -replace '[^\x20-\x7E]', '').Trim() } |
    Where-Object { $_ -match 'Ubuntu' } |
    Select-Object -First 1

if (-not $distro) {
    Write-Error "未找到 Ubuntu 发行版，请确认 WSL 中已安装 Ubuntu。"
    exit 1
}
Write-Host "检测到 WSL 发行版: $distro"

# --- 1b. 打印 verbose 列表并解析该发行版的 WSL 版本 ---
$verbose = wsl --list --verbose
Write-Host "-----------------------------------------"
Write-Host "wsl --list --verbose 输出:"
$verbose | ForEach-Object { ($_ -replace '[^\x20-\x7E]', '').Trim() } | Where-Object { $_ } | ForEach-Object { Write-Host $_ }
Write-Host "-----------------------------------------"

$wslVersion = $null
$verbose | ForEach-Object { ($_ -replace '[^\x20-\x7E]', '') -replace '\*','' } | ForEach-Object {
    $parts = ($_ -split '\s+') | Where-Object { $_ }
    if ($parts.Count -ge 3 -and $parts[0] -eq $distro) {
        $script:wslVersion = $parts[-1]
    }
}
Write-Host "该发行版 WSL 版本: $wslVersion"

if ($wslVersion -eq '1') {
    Write-Warning "该发行版是 WSL1，没有 VHDX 虚拟磁盘，无法用 compact 压缩空间。"
    Write-Host "请先将其转换为 WSL2:  wsl --set-version $distro 2"
    exit 1
}

# --- 2. 从注册表读取该发行版的 BasePath，定位 ext4.vhdx ---
$lxss    = 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Lxss'
$entries = Get-ItemProperty -Path "$lxss\*" -ErrorAction SilentlyContinue
$entry   = $entries |
    Where-Object { $_.DistributionName -eq $distro -and $_.Version -eq 2 } |
    Select-Object -First 1

if (-not $entry) {
    Write-Error "未在注册表找到 $distro 的 WSL2 配置(BasePath)。"
    exit 1
}
$vhdxPath = Join-Path $entry.BasePath 'ext4.vhdx'
if (-not (Test-Path $vhdxPath)) {
    Write-Error "注册表指向的镜像文件不存在: $vhdxPath"
    exit 1
}
$vhdx = Get-Item $vhdxPath

# --- 打印镜像位置 ---
Write-Host "========================================="
Write-Host "WSL Ubuntu 镜像位置:"
Write-Host $vhdx.FullName
Write-Host ("当前大小: {0:N2} GB" -f ($vhdx.Length / 1GB))
Write-Host "========================================="

# --- 确认门 1：是否执行「释放」(WSL 内 fstrim) ---
$decision1 = Confirm-Step "是否执行 释放(fstrim) 操作?"
if ($decision1 -eq 'n') {
    Write-Host "已取消，脚本结束。"
    exit 0
} elseif ($decision1 -eq 's') {
    Write-Host "跳过释放步骤。"
} else {
    Write-Host "步骤1: 在 WSL 内部执行 fstrim，释放已删除文件占用的块..."
    wsl -d $distro -u root -- fstrim -v /
}

# --- 确认门 2：是否执行「压缩」(写磁盘) ---
$decision2 = Confirm-Step "是否执行 VHDX 压缩 (compact vdisk)?"
if ($decision2 -eq 'n') {
    Write-Host "已取消，脚本结束。"
    exit 0
}
# 无论上一步是否跳过，压缩前都先确保 WSL 已停止以解锁镜像
Write-Host "关闭 WSL 以解锁镜像文件 (会终止所有发行版)..."
wsl --shutdown
Start-Sleep -Seconds 5

if ($decision2 -eq 's') {
    Write-Host "跳过压缩步骤。"
} else {
    Write-Host "步骤3: 压缩 VHDX 镜像..."
    if (Get-Command Optimize-VHD -ErrorAction SilentlyContinue) {
        Write-Host "使用 Optimize-VHD..."
        Optimize-VHD -Path $vhdx.FullName -Mode Full
    } else {
        Write-Host "使用 diskpart 压缩 (Optimize-VHD 不可用)..."
        $diskpart = @"
select vdisk file="$($vhdx.FullName)"
attach vdisk readonly
compact vdisk
detach vdisk
"@
        $diskpart | diskpart
    }
}

# --- 完成：打印压缩后大小 ---
$after = Get-Item $vhdx.FullName
Write-Host "========================================="
Write-Host ("压缩后大小: {0:N2} GB" -f ($after.Length / 1GB))
Write-Host "========================================="

# --- 确认门 3：是否备份 vhdx 到当前目录 ---
$decision3 = Confirm-Step "是否将 VHDX 镜像备份到当前目录?"
if ($decision3 -eq 'n') {
    Write-Host "已取消，脚本结束。"
    exit 0
} elseif ($decision3 -eq 's') {
    Write-Host "跳过备份步骤。"
} else {
    $dateStr    = Get-Date -Format 'yyyy-MM-dd'
    $backupName = "$($distro)_$dateStr.vhdx"
    $backupPath = Join-Path (Get-Location) $backupName
    Write-Host "备份中: $($vhdx.FullName)"
    Write-Host "   -> $backupPath"
    Copy-Item -Path $vhdx.FullName -Destination $backupPath -Force
    $backupSize = (Get-Item $backupPath).Length
    Write-Host ("备份完成，文件大小: {0:N2} GB" -f ($backupSize / 1GB))
    Write-Host "备份路径: $backupPath"
}

Write-Host "完成。"
