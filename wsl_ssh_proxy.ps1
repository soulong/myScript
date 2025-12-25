# ========== 日志配置 ==========
$logPath = "F:\wsl-ssh-proxy.log"

function Write-Log {
    param([string]$Message)
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    "$timestamp - $Message" | Out-File -FilePath $logPath -Encoding UTF8 -Append
}

# ========== 启动 WSL ==========
Write-Log "Try Start WSL..."
Start-Process -FilePath "wsl" -ArgumentList "-d Ubuntu-24.04", "sleep", "infinity" -WindowStyle Hidden

# ========== 等待并获取 IP ==========
$wslIp = ""
$retry = 0
$maxRetries = 20

while ($wslIp -eq "" -and $retry -lt $maxRetries) {
    Start-Sleep -Seconds 1
    $wslIp = (wsl hostname -I).Trim()
    $retry++
}

if ($wslIp -eq "") {
    Write-Log "WSL Failed or Can't Get IP"
    exit 1
}

Write-Log "WSL IP: $wslIp"

# 删除旧规则
netsh interface portproxy delete v4tov4 listenport=2222 listenaddress=0.0.0.0 > $null 2>&1

# 添加新规则
netsh interface portproxy add v4tov4 listenport=2222 listenaddress=0.0.0.0 connectport=2222 connectaddress=$wslIp

# 防火墙
$rule = Get-NetFirewallRule -DisplayName "WSL2 SSH" -ErrorAction SilentlyContinue
if (-not $rule) {
    New-NetFirewallRule -DisplayName "WSL2 SSH" -Direction Inbound -Protocol TCP -LocalPort 2222 -Action Allow | Out-Null
    Write-Log "Create Firewall Rule 'WSL2 SSH'"
}

Write-Log "Set Port Forwarding: Windows:2222 → WSL($wslIp):2222"



# # ========== 如何创立计划任务， 在PowerShell中以管理员运行 ==========
# # 创建任务：系统启动后 30 秒运行（避免 WSL 服务未就绪）
# $scriptPath = "F:\GitHub\myScript\wsl_ssh_proxy.ps1"
# $taskName = "Start WSL SSH Proxy"
# schtasks /Create /TN "$taskName" /TR "PowerShell.exe -NoProfile -ExecutionPolicy Bypass -File `"$scriptPath`"" /SC ONSTART /DELAY 0000:30 /RL HIGHEST /F