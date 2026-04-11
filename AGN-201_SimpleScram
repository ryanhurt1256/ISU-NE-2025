#
# PURPOSE
# -------
# A panel-realistic reactor simulator display:
#   • Runs 1-group point kinetics with 1 delayed-neutron precursor group.
#   • Maps a single Fine Control Rod (FCR) position → reactivity ρ(t).
#   • Computes neutron state n(t), precursor C(t), then maps n → power (W).
#   • Synthesizes two instrument channels (CH2 log, CH3 linear) from power.
#   • Allows on-screen controls, and supports physical controls via serial.
#   • Provides SCRAM behaviors (manual, external, power switch) and auto-SCRAM on 0 W / 10 W.
#
# NOTE ON MODELING CHOICES
# ------------------------
# This is a simplified dynamics model:
#   • One delayed group (instead of 6). Good for “feel” and stability, not for licensing work.
#   • Uses an instantaneous prompt-jump correction when reactivity changes.
#   • Uses RK4 integration with sub-stepping for delayed evolution stability at small Λ.
#   • Instrument channels use first-order lag for display smoothing.
#
# SIGN CONVENTION (IMPORTANT)
# ---------------------------
# ρ > 0  : supercritical (power tends to rise)
# ρ = 0  : critical (steady power if exactly balanced)
# ρ < 0  : subcritical (power tends to fall)
#
# FCR MAPPING
# -----------
# Rod withdrawn (down) → more negative ρ
# Rod inserted (up)  → more positive ρ
#
# TRANSPORT INPUT FORMAT
# ------------------------------
# This build supports either USB serial or Ethernet using the SAME Arduino packet format:
#   • telemetry packet: packet1, packet2, posSetH, posSetL, posReadH, posReadL, knob1H, knob1L, knob2H, knob2L, 0x24
#   • optional ASCII debug lines may appear on the same stream and are skipped when the payload is implausible
#
#     rod_frac: positionRead / 1023.0   (0 withdrawn, 1 inserted)
#     ch2_det : quantized knob1 → -3..+3
#     ch3_det : quantized knob2 → -3..+3
#     power_sw: bit 1 of packet1 (latched)
#     scram_sw: bit 0 of packet1 (latched)
#
# Outbound command to the Arduino sketch over either transport:
#   bytes([0x24, 0xFF, flipMask])
#     bit 0 -> toggle SCRAM latch
#     bit 1 -> toggle POWER latch
#
# Assumption for Ethernet mode:
#   the Arduino Ethernet firmware sends and receives the same byte stream as the USB serial firmware.
# ------------------------------------------------------------------------------

import pygame                     # Pygame for windowing, drawing, keyboard/mouse input, timing
import sys                        # sys.exit() for clean termination
import math                       # math for trig/log10/hypot conversions and geometry
import time                       # time.time() for wall-clock dt and scram timing
import collections                # collections.deque for fixed-length plot histories
import socket                     # Bluetooth RFCOMM (disabled by default in this Arduino-matched build)

# ---------- OPTIONAL SERIAL ----------
try:
    import serial                 # pyserial: reads hardware controls from Arduino/USB serial
except Exception:
    serial = None                 # if pyserial not installed, disable serial I/O cleanly

# =========================
# PHYSICS CONFIG (TUNED 3/14/2026)
# =========================
beta   = 0.00745                  # βeff: effective delayed neutron fraction from AGN-201 appendix (dimensionless)
Lambda = 6.22e-5                  # Λ: neutron generation lifetime from AGN-201 appendix (s)

lam    = 0.3744                   # λeff: 1-group delayed neutron decay constant (1/s), obtained by averaging 6-group data to one group using the β-weighted average λeff = Σ(βi*λi)/Σβi

# Prompt jump model
PROMPT_JUMP_ENABLE = True         # apply 1-group prompt-jump correction when applied reactivity changes
PROMPT_JUMP_EPS    = 1e-8         # minimum |Δρ| required before applying prompt jump (prevents tiny numerical chatter)

# Fully withdrawn is negative; inserting adds + worth
RHO_WITHDRAWN   = -7.44e-4        # ρ at fully withdrawn rod position (subcritical baseline), chosen so that ~10 cm inserted is near critical, consistent with lab procedure
RHO_WORTH_MAG   = +3.10e-3        # total positive worth of the fuel-loaded FCR, 0.00310 Δk/k

# Hard safety bounds
RHO_NEG_LIMIT = -0.20             # clamp lower bound to avoid runaway negative stiffness / nonsense
RHO_POS_LIMIT = 0.90 * beta       # clamp upper bound below β to avoid prompt-supercritical blowup

# Reactivity smoothing
RHO_TAU      = 0.00               # time constant (s) for lag filter on commanded ρ → applied ρ
RHO_DEADBAND = 2e-5               # small deadband to prevent chatter when changes are tiny

# Power mapping
MAX_POWER_DISPLAY = 5.0           # plot-only cap for graph readability (does NOT limit physics)

# =========================
# AUTO-SCRAM
# =========================
SCRAM_HIGH_W = 10.0               # auto-scram if computed power >= 10 W
SCRAM_LOW_W  = 1e-6               # “0 W” threshold: scram if power stays at numerical floor
SCRAM_DEBOUNCE_S = 0.25           # time (s) power must violate limit before SCRAM triggers

# =========================
# INSTRUMENT CONFIG
# =========================
I_FLOOR = 1e-13                   # minimum current to avoid log10(0) and represent leakage floor

CH3_HIGH_SCRAM_A = 1.8e-7         # reference: current at 10 W full-scale on CH3 high range (A)
I_PER_WATT = CH3_HIGH_SCRAM_A / 10.0  # linear conversion (A/W) based on reference above
I_PER_WATT_CH2 = I_PER_WATT       # channel 2 current-per-watt (same scaling here)
I_PER_WATT_CH3 = I_PER_WATT       # channel 3 current-per-watt (same scaling here)

CH2_TAU = 0.25                    # instrument lag (s) for CH2 indication smoothing
CH3_TAU = 0.25                    # instrument lag (s) for CH3 indication smoothing

CH2_LOG_SPAN_DECADES = 4          # CH2 shows a 4-decade log window at any dial setting
CH2_LOG_MAX_AT_0     = -8         # at detent=0, top of CH2 window is 10^-8 A (then shifts by detent)

CH3_RANGES = [                    # CH3 full-scale current ranges (A), detent -3..+3 maps to these
    0.3e-12, 0.3e-11, 0.3e-10,
    0.3e-9,  0.3e-8,  0.3e-7,
    0.3e-6
]

# =========================
# FCR SHAPE CURVE
# =========================
FCR_TRAVEL_CM = 23.0              # physical rod travel (cm) used to convert fraction ↔ cm

FCR_CURVE_POINTS_CM_PCT = [       # (cm_inserted, %total_worth_inserted) points from FCR worth curve
    (0.0,   0.0),
    (1.0,   0.5),
    (2.0,   1.5),
    (3.0,   2.7),
    (4.0,   4.0),
    (5.0,   5.5),
    (6.0,   8.0),
    (7.0,  11.0),
    (8.0,  15.0),
    (9.0,  19.0),
    (10.0, 24.0),
    (11.0, 30.0),
    (12.0, 36.0),
    (13.0, 43.0),
    (14.0, 50.0),
    (15.0, 57.0),
    (16.0, 65.0),
    (17.0, 72.0),
    (18.0, 80.0),
    (19.0, 87.0),
    (20.0, 94.0),
    (21.0, 99.0),
    (22.0, 100.0),
    (23.0, 100.0),
]

# =========================
# UI CONFIG
# =========================
pygame.init()                     # initialize all Pygame subsystems (display, font, input, etc.)

#WIDTH, HEIGHT = 1700, 1200        # window pixel size
#screen = pygame.display.set_mode((WIDTH, HEIGHT))  # create window surface
screen = pygame.display.set_mode((0, 0), pygame.FULLSCREEN) # Use native display resolution (important for Raspberry Pi)
WIDTH, HEIGHT = screen.get_size() # window pixel size of display being used
pygame.display.set_caption("Control Panel (FCR) — CH2 LOG + CH3 LIN")  # window title
clock = pygame.time.Clock()       # FPS timing controller (used with clock.tick)

WHITE  = (255, 255, 255)          # RGB color constants for drawing
GREY   = (150, 150, 150)
BLUE   = (50, 100, 255)
RED    = (255, 60, 60)
YELLOW = (255, 255, 0)
BLACK  = (0, 0, 0)
GREEN  = (0, 200, 0)
SCRAM_RED = (255, 0, 0)
DIAL_INACTIVE_COLOR  = (50, 50, 50)    # dial face color when idle
LEVER_INACTIVE_COLOR = (80, 80, 80)    # lever handle color when inactive / disabled

# ---------- FONT CACHE ----------
# Create fonts once at startup
font = pygame.font.SysFont(None, 24)       # medium font for key header text
inst_font = pygame.font.SysFont(None, 16)  # smaller font for instrument/status lines
title_font = pygame.font.SysFont(None, 18) # small title font for plot labels
tab_font = pygame.font.SysFont(None, 22)   # font for collapse/expand tab
dial_tick_font = pygame.font.SysFont(None, 18)  # small font for dial tick labels

# ---------- BASE LAYOUT ----------
MARGIN_LEFT, MARGIN_RIGHT, MARGIN_BOTTOM = 30, 30, 20  # margins; bottom mostly reserved for toggle tab

GRAPH_TOP  = 10                 # top padding above header/graph stack
HEADER_H   = 100                # vertical space reserved for header text
GRAPH_X    = MARGIN_LEFT        # left edge of plot region
GRAPH_Y    = GRAPH_TOP + HEADER_H   # top of plot region (below header)
GRAPH_WIDTH  = WIDTH - MARGIN_LEFT - MARGIN_RIGHT  # usable horizontal plot width

BASE_GRAPH_HEIGHT = int(HEIGHT * 0.34)  # default plot height when controls are visible
TOGGLE_H = 30                   # height of collapse/expand “Controls” tab
TOGGLE_GAP = 8                  # gap between graph region and the toggle tab

# Controls collapse state
controls_collapsed = False      # False = controls visible; True = controls hidden and graphs expand

# Dynamic layout vars
GRAPH_HEIGHT = BASE_GRAPH_HEIGHT  # actual current graph height (changes if collapsed)
linear_width = GRAPH_WIDTH // 2   # each plot gets half the width
right_start  = GRAPH_X + linear_width  # x origin for the right-hand plot

# Lever geometry
LEVER_TRACK_TOP    = 0           # y (px) of lever track top; computed in update_layout
LEVER_TRACK_BOTTOM = 0           # y (px) of lever track bottom; computed in update_layout
LEVER_HEIGHT       = 50          # lever handle height (px)
LEVER_X            = GRAPH_X + GRAPH_WIDTH // 2  # lever centered under plots
lever_y            = 0           # current lever handle y-position (px); updated later
lever_name = "Fine CR"           # label shown under lever

# Dials
NUM_DIALS   = 2                  # number of dials (CH2 and CH3)
DIAL_RADIUS = 38                 # dial radius (px)
dial_centers = [(0, 0), (0, 0)]  # dial center positions; computed in update_layout
dial_names  = ["CH2 Range (LOG)", "CH3 Range (LIN)"]  # labels
dial_detents = [0, 0]            # dial detents (integers -3..+3), applied to CH2 and CH3

# Buttons
button_width, button_height = 100, 50   # size for POWER button
power_on = False                        # system “powered” state; gates physics updates
#POWER_RECT = pygame.Rect(60, HEIGHT - 155, button_width, button_height)  # POWER button rectangle
POWER_RECT = pygame.Rect(0, 0, button_width, button_height)  # placeholder; real position set in update_layout()
SCRAM_RECT = pygame.Rect(70, 0, 90, 90)  # SCRAM button rectangle; y fixed in update_layout

# Toggle tab rect
TOGGLE_RECT = pygame.Rect(GRAPH_X, GRAPH_Y + GRAPH_HEIGHT + TOGGLE_GAP, GRAPH_WIDTH, TOGGLE_H)  # placeholder, updated later

# Selection + cooldown
selected_lever = False            # whether the lever has been selected for keyboard control
selected_dial  = None             # index of selected dial (0 or 1), or None
lever_active = False              # lever “active” (lit) state; times out after inactivity
dial_active  = [False, False]     # per-dial “active” state; times out after inactivity
lever_activation_time = 0.0       # last time lever interaction occurred
dial_activation_time  = [0.0, 0.0]# last time each dial interaction occurred
COOLDOWN_TIME = 3.0               # seconds to auto-deselect controls after last use

scram_active = False              # SCRAM in progress flag
scram_start_time = 0.0            # when SCRAM started (wall-clock)
scram_reason = ""                 # readable reason displayed on-screen

# Histories (for plots)
ch2_history = collections.deque(maxlen=300)   # CH2 time history buffer (most recent 300 samples)
ch3_history = collections.deque(maxlen=300)   # CH3 time history buffer
power_history = collections.deque(maxlen=300) # power time history buffer (currently not drawn, but retained)

# Indicated currents
ch2_I_ind = I_FLOOR               # current displayed by CH2 after lag filtering
ch3_I_ind = I_FLOOR               # current displayed by CH3 after lag filtering

# Kinetics state
n = 1e-12                         # neutron/power state amplitude; here displayed power is taken directly from n
C = 0.0                           # delayed neutron precursor concentration (arbitrary units)
rho_prev = 0.0                    # previous applied reactivity (for filtering and continuity)
rho_now  = 0.0                    # current applied reactivity used in kinetics RHS

# Time / neutron count display
sim_time_s = 0.0                  # simulation time accumulator (only increments when power_on and not scram)
power_on_wall_time = 0.0          # wall-clock timestamp when power was turned on (informational / future use)

# Auto-scram debounce timers
high_scram_timer = 0.0            # accumulates time spent above SCRAM_HIGH_W
low_scram_timer  = 0.0            # accumulates time spent below SCRAM_LOW_W

last_time = time.time()           # previous frame wall-clock time used to compute dt

# ---------- TRANSPORT (USB serial or Ethernet) ----------
BAUDRATE = 9600                   # serial speed; must match Arduino sketch
SERIAL_PORTS = (                  # candidate device paths to try (Linux/RPi naming)
    "/dev/ttyACM0", "/dev/ttyACM1",
    "/dev/ttyUSB0", "/dev/ttyUSB1",
    "/dev/serial0"
)
TRANSPORT_MODE = "auto"          # "auto", "serial", or "ethernet"
TRANSPORT_RETRY_S = 2.0           # retry interval if the Arduino link is not ready when Python starts
ETH_HOST = "192.168.1.50"        # change to the Arduino Ethernet IP when Ethernet firmware is enabled
ETH_PORT = 5000                   # change to the Arduino Ethernet TCP port when Ethernet firmware is enabled
ETH_CONNECT_TIMEOUT_S = 0.20      # short timeout keeps GUI responsive during retries

transport = None                  # dict describing the active link: {"kind":..., "obj":..., "name":...}
transport_name = "DISCONNECTED"  # human-readable link status for the header
transport_next_retry = 0.0        # next wall-time when a reconnect attempt is allowed
transport_rx_packets = 0          # good telemetry packets parsed from the Arduino
transport_tx_commands = 0         # outbound GUI command packets sent to the Arduino
transport_last_rx_time = 0.0      # wall-time of the most recent valid telemetry packet
transport_last_tx_time = 0.0      # wall-time of the most recent outbound command
transport_last_error = ""        # last transport open/write/read error string for debugging


def setup_serial_transport():
    global transport_name, transport_last_error
    if serial is None:
        transport_name = "PYSerial missing"
        return None
    for p in SERIAL_PORTS:
        try:
            s = serial.Serial(p, BAUDRATE, timeout=0, write_timeout=0)
            time.sleep(2.0)       # opening USB serial resets many Arduinos; let the bootloader finish
            try:
                s.reset_input_buffer()
                s.reset_output_buffer()
            except Exception:
                pass
            transport_name = p
            transport_last_error = ""
            return {"kind": "serial", "obj": s, "name": p}
        except Exception as exc:
            transport_last_error = str(exc)
    transport_name = "DISCONNECTED"
    return None


def setup_ethernet_transport():
    global transport_name, transport_last_error  # allow function to update human-readable link status + last error text
    s = None                                      # predefine socket variable so it can be safely closed in the exception path

    try:
        s = socket.socket(
            socket.AF_INET,                       # IPv4 Internet socket
            socket.SOCK_STREAM                    # TCP stream socket (reliable ordered byte stream)
        )
        s.settimeout(ETH_CONNECT_TIMEOUT_S)       # short connect timeout prevents the GUI from freezing during retries
        s.setsockopt(
            socket.IPPROTO_TCP,                   # apply option at the TCP protocol layer
            socket.TCP_NODELAY,                   # disable Nagle's algorithm so small control packets send immediately
            1
        )
        s.connect((ETH_HOST, ETH_PORT))          # connect outward to the Arduino Ethernet server at the configured IP/port
        s.setblocking(False)                     # switch to non-blocking mode so later recv() calls do not stall the frame loop

        name = f"tcp://{ETH_HOST}:{ETH_PORT}"    # readable transport name shown in the simulator header/status line
        transport_name = name                    # store connected transport name for display/debugging
        transport_last_error = ""                # clear any previous link error because connection succeeded

        return {
            "kind": "ethernet",                  # identifies this transport as Ethernet/TCP
            "obj": s,                            # actual connected socket object used for send/recv
            "name": name                         # readable name for status text
        }

    except Exception as exc:
        transport_last_error = str(exc)          # save the exception text so the UI can show the most recent connection failure

        if s is not None:                        # only try to close if the socket object was successfully created
            try:
                s.close()                        # clean up partially opened socket so retries start fresh
            except Exception:
                pass                             # ignore close errors because connection already failed and we are recovering

        transport_name = "DISCONNECTED"          # update visible status to show that no Ethernet link is active
        return None                              # tells the caller that Ethernet setup failed this attempt


def open_transport():
    # Auto mode prefers Ethernet first when available, then falls back to USB serial.
    if TRANSPORT_MODE == "ethernet":
        return setup_ethernet_transport()
    if TRANSPORT_MODE == "serial":
        return setup_serial_transport()
    tr = setup_ethernet_transport()
    if tr is not None:
        return tr
    return setup_serial_transport()


def ensure_transport_connected(now):
    global transport, transport_next_retry
    if transport is not None:
        return
    if now < transport_next_retry:
        return
    transport_next_retry = now + TRANSPORT_RETRY_S
    transport = open_transport()


def close_transport():
    global transport, transport_name
    if transport is not None:
        try:
            transport["obj"].close()
        except Exception:
            pass
    transport = None
    transport_name = "DISCONNECTED"


def write_transport_bytes(transport_obj, payload):
    if transport_obj is None:
        return False
    if transport_obj["kind"] == "serial":
        transport_obj["obj"].write(payload)
        transport_obj["obj"].flush()
        return True
    if transport_obj["kind"] == "ethernet":
        transport_obj["obj"].sendall(payload)
        return True
    raise RuntimeError("Unknown transport kind")


def read_transport_bytes(transport_obj):
    if transport_obj is None:
        return b""
    if transport_obj["kind"] == "serial":
        ser_obj = transport_obj["obj"]
        if hasattr(ser_obj, "in_waiting") and ser_obj.in_waiting <= 0:
            return b""
        return ser_obj.read(ser_obj.in_waiting)
    if transport_obj["kind"] == "ethernet":
        try:
            data = transport_obj["obj"].recv(4096)
            if data == b"":
                raise ConnectionError("Ethernet peer closed the socket")
            return data
        except BlockingIOError:
            return b""
    return b""


def send_flipmask_command(transport_obj, flip_mask):
    global transport_tx_commands, transport_last_tx_time, transport_last_error
    if transport_obj is None:
        return False
    try:
        write_transport_bytes(transport_obj, bytes((0x24, 0xFF, flip_mask & 0xFF)))
        transport_tx_commands += 1
        transport_last_tx_time = time.time()
        transport_last_error = ""
        return True
    except Exception as exc:
        transport_last_error = str(exc)
        close_transport()
        return False

# Persistent buffer for packet parsing (works for USB serial or Ethernet byte streams)
rx_buffer = bytearray()

# =========================
# BLUETOOTH SEND (disabled in this Arduino-tailored build)
# =========================
# The active Arduino control link is handled by the transport layer above. Bluetooth is left here but disabled
# so it cannot interfere with bench testing or be confused with the primary USB/Ethernet link.

BT_ENABLE = False               # keep False for the current Arduino-matched hardware test build
BT_ADDR = "AA:BB:CC:DD:EE:FF"   # Bluetooth MAC address of the receiving Raspberry Pi NEEDS TO BE ADJUSTED FOR BLUETOOTH DEVICE
BT_CH = 1                       # RFCOMM channel number (must match receiver)
BT_HZ = 20                      # Data transmission rate (messages per second)

bt_sock = None                  # Bluetooth socket object (None until connection established)
bt_last_send = 0.0              # Timestamp of last transmitted message
bt_next_retry = 0.0             # next time to attempt a Bluetooth connect
BT_RETRY_S = 1.5                # wait this long between connection attempts (prevents frame-freeze spam)
BT_STARTUP_DELAY_S = 1.0        # wait this long after program start before first BT connect attempt
bt_enable_time = time.time() + BT_STARTUP_DELAY_S  # earliest wall-time when BT is allowed to connect

def bt_connect():
    """Create and connect a Bluetooth RFCOMM socket (short timeout so it never stalls the GUI)."""
    s = socket.socket(                   # Create Bluetooth socket
        socket.AF_BLUETOOTH,             # Bluetooth protocol family
        socket.SOCK_STREAM,              # Stream socket (serial-like)
        socket.BTPROTO_RFCOMM            # RFCOMM protocol
    )
    s.settimeout(0.10)                   # shorter timeout further reduces visible connect hitches in the GUI loop
    s.connect((BT_ADDR, BT_CH))          # Connect to receiver Pi (can raise timeout/OS errors)
    s.settimeout(0.0)                    # set non-blocking mode for send() so Bluetooth cannot stall the frame loop
    return s                             # Return connected socket

def read_controls_from_transport(transport_obj):
    global transport_rx_packets, transport_last_rx_time, transport_last_error
    if transport_obj is None:     # if no active link, no data
        return None
    try:
        data = read_transport_bytes(transport_obj)
        if not data:
            return None

        rx_buffer.extend(data)

        # Keep the buffer bounded even if the stream gets noisy
        if len(rx_buffer) > 4096:
            del rx_buffer[:-256]

        # Process any complete packet(s) ending with 0x24. The provided Arduino sketch mixes these binary packets
        # with ASCII Serial.println() debug text on the same port, so each candidate terminator is validated before use.
        while len(rx_buffer) >= 11:
            try:
                idx = rx_buffer.index(0x24)
            except ValueError:
                break  # no terminator in the buffer yet

            if idx < 10:
                del rx_buffer[:idx + 1]      # not enough bytes before 0x24 to form a packet
                continue

            pkt = rx_buffer[idx - 10:idx]    # 10 data bytes immediately before the 0x24 terminator
            packet1 = pkt[0]
            packet2 = pkt[1]
            pos_set  = (pkt[2] << 8) | pkt[3]
            pos_read = (pkt[4] << 8) | pkt[5]
            knob1    = (pkt[6] << 8) | pkt[7]
            knob2    = (pkt[8] << 8) | pkt[9]

            # Plausibility checks tuned to the current Arduino sketch:
            #   • the four analog channels are 10-bit ADC values (0..1023)
            #   • packet2 currently only uses bit 0 (controlPin debug)
            plausible = (
                pos_set  <= 1023 and
                pos_read <= 1023 and
                knob1    <= 1023 and
                knob2    <= 1023 and
                (packet2 & 0xFE) == 0
            )

            if not plausible:
                del rx_buffer[:idx + 1]      # reject false 0x24 hits from the mixed stream and keep searching
                continue

            # === MAPPINGS FOR SIMULATOR ===
            rod_frac = pos_read / 1023.0                                 # 0..1 (physical rod feedback)
            # If your pot is wired backwards, uncomment the next line:
            # rod_frac = 1.0 - rod_frac

            # Quantize rotary knobs to detents -3..+3 (7 positions)
            ch2_det = int(round((knob1 / 1023.0) * 6 - 3))
            ch2_det = max(-3, min(3, ch2_det))
            ch3_det = int(round((knob2 / 1023.0) * 6 - 3))
            ch3_det = max(-3, min(3, ch3_det))

            power_sw = bool(packet1 & (1 << 1))   # latched power toggle bit from the Arduino sketch
            scram_sw = bool(packet1 & (1 << 0))   # latched scram toggle bit from the Arduino sketch

            del rx_buffer[:idx + 1]               # remove processed packet plus terminator from the buffer

            transport_rx_packets += 1
            transport_last_rx_time = time.time()
            transport_last_error = ""
            return (rod_frac, ch2_det, ch3_det, power_sw, scram_sw)

        return None  # no complete validated packet this frame
    except Exception as exc:
        transport_last_error = str(exc)
        close_transport()
        return None                # if parsing fails badly, drop the link and allow reconnect

# =========================
# HELPERS
# =========================
def clamp(v, lo, hi):
    return lo if v < lo else hi if v > hi else v  # saturating clamp used everywhere for safety

def first_order_lag(y, x, dt, tau):
    if tau <= 1e-6:               # if tau is ~0, treat as instantaneous
        return x
    a = dt / (tau + dt)           # discrete-time equivalent of 1st-order low-pass coefficient
    return y + a * (x - y)        # filtered output approaching x with time constant tau

def format_watts(P):
    P = float(P)                  # enforce numeric float
    if P >= 1.0:                  # format selection for readability
        return f"{P:.2f} W"
    if P >= 0.01:
        return f"{P:.4f} W"
    return f"{P:.6f} W"

def kinetics_rhs(n_s, C_s, rho):
    # POINT KINETICS (1 delayed group)
    # dn/dt = ((ρ - β)/Λ) n + λ C
    # dC/dt = (β/Λ) n - λ C
    #
    # Interpretation:
    #   • n is the amplitude of neutron population / flux / power.
    #   • C is precursor concentration producing delayed neutrons at rate λC.
    #   • ρ introduces net multiplication: ρ>0 increases dn/dt, ρ<0 decreases dn/dt.
    #   • β/Λ term moves some production into precursor inventory.
    dn = ((rho - beta) / Lambda) * n_s + lam * C_s  # neutron balance ODE
    dC = (beta / Lambda) * n_s - lam * C_s          # precursor balance ODE
    return dn, dC                                   # return derivatives for integrator

def rk4_step(n0, C0, dt, rho):
    # 4th-order Runge-Kutta (RK4) integration for stiff kinetics (Λ small).
    # RK4 reduces numerical error compared to Euler and is stable for small dt.
    k1n, k1C = kinetics_rhs(n0, C0, rho)                                        # slope at start
    k2n, k2C = kinetics_rhs(n0 + 0.5*dt*k1n, C0 + 0.5*dt*k1C, rho)             # slope at midpoint (1)
    k3n, k3C = kinetics_rhs(n0 + 0.5*dt*k2n, C0 + 0.5*dt*k2C, rho)             # slope at midpoint (2)
    k4n, k4C = kinetics_rhs(n0 + dt*k3n,     C0 + dt*k3C,     rho)             # slope at end
    n1 = n0 + (dt/6.0)*(k1n + 2*k2n + 2*k3n + k4n)                             # weighted average update
    C1 = C0 + (dt/6.0)*(k1C + 2*k2C + 2*k3C + k4C)                             # weighted average update
    return n1, C1                                                               # next state

def apply_prompt_jump(n_old, rho_old, rho_new):
    # 1-group prompt-jump approximation:
    #   n+ / n- = (β - ρ_old) / (β - ρ_new)
    #
    # Assumptions:
    #   • precursor concentration C does not change instantaneously
    #   • reactivity change is effectively instantaneous to the neutron population
    #   • valid only while ρ remains below prompt critical (ρ < β)
    denom = beta - rho_new                        # denominator of prompt-jump ratio
    if denom <= 1e-12:                            # extra safety guard against divide-by-zero / nonsense
        return n_old                              # if somehow invalid, leave neutron state unchanged

    numer = beta - rho_old                        # numerator of prompt-jump ratio
    if numer <= 1e-12:                            # extra guard against invalid old state near/above prompt critical
        return n_old

    ratio = numer / denom                         # prompt-jump multiplier from old to new reactivity
    ratio = clamp(ratio, 0.1, 10.0)              # keep single-frame jump bounded against numerical spikes
    return max(1e-12, n_old * ratio)             # apply jump and keep neutron state above numerical floor

def smooth_points(points, steps=6): # Simple linear interpolation densification used to make CH3 line look smoother
    if len(points) < 2:             # not enough points to smooth
        return points
    out = []                        # output list of densified points
    for i in range(len(points) - 1):# for each segment
        x1, y1 = points[i]          # start of segment
        x2, y2 = points[i + 1]      # end of segment
        for t in range(steps):      # insert intermediate points
            a = t / steps           # interpolation fraction
            out.append((x1*(1-a) + x2*a, y1*(1-a) + y2*a))  # linear interpolation
    out.append(points[-1])          # include final original point
    return out

def angle_from_detent(det):
    return -90.0 + 30.0 * float(clamp(det, -3, 3))  # detent -3..3 maps to pointer angle degrees

def draw_dial_ticks(cx, cy, radius):
    angle_step = 30                                 # degrees between detents
    for i in range(7):                              # 7 detents: -3,-2,-1,0,+1,+2,+3
        det = i - 3                                 # convert loop index to detent value
        angle_deg = -90 - angle_step * det          # tick angle
        ang = math.radians(angle_deg)               # convert degrees to radians for trig
        tick_len = 10 if det == 0 else 6            # longer tick at center detent for emphasis

        ox = cx + radius * math.cos(ang)            # outer tick endpoint x
        oy = cy + radius * math.sin(ang)            # outer tick endpoint y
        ix = cx + (radius - tick_len) * math.cos(ang)  # inner tick endpoint x
        iy = cy + (radius - tick_len) * math.sin(ang)  # inner tick endpoint y
        pygame.draw.line(screen, BLACK, (ox, oy), (ix, iy), 2)  # draw tick mark

        lx = cx + (radius + 22) * math.cos(ang)     # label position x (outside dial)
        ly = cy + (radius + 22) * math.sin(ang)     # label position y
        label_txt = "0" if det == 0 else str(det)   # print 0 at center, otherwise detent number
        lab = dial_tick_font.render(label_txt, True, BLACK)   # rasterize label text surface
        screen.blit(lab, lab.get_rect(center=(lx, ly)))  # blit centered at label position

def lever_inserted_fraction(y):
    travel = (LEVER_TRACK_BOTTOM - LEVER_TRACK_TOP - LEVER_HEIGHT)  # total usable travel (px)
    if travel <= 0:                        # safety: avoid divide-by-zero if layout invalid
        return 0.0
    x = (LEVER_TRACK_BOTTOM - LEVER_HEIGHT - y) / travel  # map y to fraction (top=1 inserted, bottom=0 withdrawn)
    return clamp(x, 0.0, 1.0)            # clamp to 0..1 physical fraction

def lever_y_from_inserted(x):
    travel = (LEVER_TRACK_BOTTOM - LEVER_TRACK_TOP - LEVER_HEIGHT)  # total usable travel (px)
    return (LEVER_TRACK_BOTTOM - LEVER_HEIGHT) - clamp(x, 0.0, 1.0) * travel  # fraction → y position

def interp_piecewise(points, x):
    if not points:                         # if no data points, return 0
        return 0.0
    if x <= points[0][0]:                  # below first x, clamp to first y
        return points[0][1]
    if x >= points[-1][0]:                 # above last x, clamp to last y
        return points[-1][1]
    for i in range(len(points) - 1):       # search which segment contains x
        x0, y0 = points[i]                 # segment start
        x1, y1 = points[i + 1]             # segment end
        if x0 <= x <= x1:                  # if x is inside this interval
            t = (x - x0) / max(1e-12, (x1 - x0))  # normalized fraction (avoid division by tiny)
            return y0 + t * (y1 - y0)      # linear interpolate y
    return points[-1][1]                   # fallback (should not hit)

def fcr_inserted_worth_fraction(x_inserted):                # Convert inserted fraction x_inserted ∈ [0,1] to “fraction of total worth inserted”.  Table is “% total worth vs cm inserted”.
    x_inserted = clamp(float(x_inserted), 0.0, 1.0)         # enforce numeric and clamp
    cm_inserted = x_inserted * FCR_TRAVEL_CM                # convert inserted fraction to inserted depth (cm)
    pct = interp_piecewise(FCR_CURVE_POINTS_CM_PCT, cm_inserted)  # interpolate % inserted-worth from table
    f_inserted = clamp(pct / 100.0, 0.0, 1.0)               # convert to fraction 0..1
    return f_inserted                                       # final fraction of total worth inserted

def reset_on_power_on(now): # Reset state variables when transitioning from OFF → ON (fresh run condition).
    global n, C, rho_prev, rho_now
    global ch2_I_ind, ch3_I_ind
    global sim_time_s, power_on_wall_time
    global high_scram_timer, low_scram_timer  

    power_history.clear()          # clear stored power trace
    ch2_history.clear()            # clear CH2 trace
    ch3_history.clear()            # clear CH3 trace

    sim_time_s = 0.0               # reset simulation time counter
    power_on_wall_time = now       # remember wall time when power was enabled

    n = 1.0                        # initialize state so displayed power starts at about 1 W (here power = n)
    C = (beta / (Lambda * lam)) * n# set C near the steady-state relation for 1-group kinetics at ρ≈0
    rho_prev = 0.0                 # reset applied reactivity
    rho_now = 0.0                  # reset applied reactivity

    P0 = 1.0                       # “initial power” for instrument initialization (W)
    ch2_I_ind = max(I_FLOOR, I_PER_WATT_CH2 * P0)  # initialize CH2 indicated current
    ch3_I_ind = max(I_FLOOR, I_PER_WATT_CH3 * P0)  # initialize CH3 indicated current

    high_scram_timer = 0.0         # reset high-power debounce timer
    low_scram_timer  = 0.0         # reset low-power debounce timer

def start_scram(reason, now):
    global scram_active, scram_start_time, scram_reason
    global high_scram_timer, low_scram_timer
    if not scram_active:           # only start SCRAM once per event
        scram_active = True        # latch SCRAM active
        scram_start_time = now     # mark start time
        scram_reason = reason      # store reason for display
    high_scram_timer = 0.0         # clear debounce timers so old state cannot retrigger
    low_scram_timer  = 0.0

def update_layout(collapsed):
    """
    layout:
      • Graph grows when controls are collapsed.
      • When controls are visible, controls are packed into the remaining screen height.
      • Keeps lever, SCRAM, POWER, and dials on-screen for smaller Raspberry Pi displays.
    """
    global GRAPH_HEIGHT, linear_width, right_start
    global LEVER_TRACK_TOP, LEVER_TRACK_BOTTOM, lever_y
    global SCRAM_RECT, POWER_RECT, dial_centers, TOGGLE_RECT

    xins = lever_inserted_fraction(lever_y) if LEVER_TRACK_BOTTOM > LEVER_TRACK_TOP else 0.0  # preserve current lever inserted fraction so layout does not change rod position 

    if collapsed:
        bottom_reserve = TOGGLE_GAP + TOGGLE_H + 12  # reserve space below the enlarged graph for the always-visible hide/show tab
        GRAPH_HEIGHT = max(220, HEIGHT - (GRAPH_Y + bottom_reserve))  # expand graph downward when controls are hidden, but keep a minimum height
    else:
        GRAPH_HEIGHT = min(int(HEIGHT * 0.34), 260)  # Graph height

    linear_width = GRAPH_WIDTH // 2  # split total graph width evenly between CH2 and CH3
    right_start  = GRAPH_X + linear_width  # x-position where the right-hand plot begins

    TOGGLE_RECT = pygame.Rect(
        GRAPH_X,                                  # align tab with graph left edge
        GRAPH_Y + GRAPH_HEIGHT + TOGGLE_GAP,      # place tab directly below the graphs
        GRAPH_WIDTH,                              # tab spans full graph width
        TOGGLE_H                                  # tab height
    )

    if collapsed:  # Keep control geometry valid even while hidden so helper functions still have safe coordinates
        LEVER_TRACK_TOP = TOGGLE_RECT.bottom + 10
        LEVER_TRACK_BOTTOM = LEVER_TRACK_TOP + 120
        lever_y = lever_y_from_inserted(xins)  # restore lever position using the preserved inserted fraction
        SCRAM_RECT = pygame.Rect(70, LEVER_TRACK_TOP + 10, 90, 90)  # temporary off-focus SCRAM placement while controls are collapsed
        POWER_RECT = pygame.Rect(20, HEIGHT - 75, button_width, button_height)  # temporary POWER placement while controls are collapsed
        dial_centers = [(0, 0), (0, 0)]  # hide dial centers while controls are collapsed
        return

    # -------- Controls region below toggle tab --------
    controls_top = TOGGLE_RECT.bottom + 10   # start controls just below the tab with a small visual gap
    controls_bottom = HEIGHT - 10            # keep a small bottom margin so controls do not touch the screen edge
    controls_h = max(140, controls_bottom - controls_top)  # total vertical space available for the visible controls block

    track_h = min(120, max(80, controls_h - 70))  # compact lever travel height chosen to balance visibility and available space
    LEVER_TRACK_TOP = controls_top + 5            # top of lever track within the controls band
    LEVER_TRACK_BOTTOM = LEVER_TRACK_TOP + track_h  # bottom of lever track based on compact Pi-friendly track height

    lever_y = lever_y_from_inserted(xins)  # restore lever handle position after layout changes without changing rod insertion

    # SCRAM button placed left of lever region
    SCRAM_RECT = pygame.Rect(25, LEVER_TRACK_TOP + max(0, (track_h - 90) // 2), 90, 90)

    # POWER button near bottom-left
    POWER_RECT = pygame.Rect(20, SCRAM_RECT.bottom + 18, button_width, button_height)

    # Dials placed below the lever region while still staying on-screen
    dial_y = min(HEIGHT - 60, LEVER_TRACK_BOTTOM + 45)  # final dial y-position limited so the dials never fall below the display
    dial_spacing = GRAPH_WIDTH // (NUM_DIALS + 1)       # spread the two dials evenly across the available control width
    dial_centers = [
        (GRAPH_X + dial_spacing * (i + 1), dial_y)      # compute each dial center from equal spacing and the shared dial y-position
        for i in range(NUM_DIALS)
    ]
lever_y = GRAPH_Y + BASE_GRAPH_HEIGHT                 # temporary y so lever_inserted_fraction won't divide by zero
LEVER_TRACK_TOP = GRAPH_Y + BASE_GRAPH_HEIGHT + 95    # set initial lever track top
LEVER_TRACK_BOTTOM = LEVER_TRACK_TOP + 270            # set initial lever track bottom
lever_y = LEVER_TRACK_BOTTOM - LEVER_HEIGHT           # start lever fully withdrawn (down)
update_layout(controls_collapsed)                     # compute final positions from layout state

# =========================
# MAIN LOOP
# =========================
running = True                                        # main loop flag
while running:
    screen.fill(WHITE)                                # clear entire screen each frame

    now = time.time()                                 # current wall-clock time
    ensure_transport_connected(now)                   # keep retrying until the Arduino transport appears
    dt = now - last_time                              # frame time step
    if dt <= 0.0:                                     # guard against zero/negative dt
        dt = 1.0/60.0                                 # assume ~60 FPS step
    if dt > 0.1:                                      # clamp dt to avoid huge physics jumps on pauses
        dt = 0.1
    last_time = now                                   # update last_time for next frame

    # ---------- SERIAL INPUT ----------
    sc = read_controls_from_transport(transport)      # attempt to read one telemetry packet
    if sc is not None:                                # if valid data received
        rod_frac, ch2_det, ch3_det, power_sw, scram_sw = sc  # unpack controls
        lever_y = lever_y_from_inserted(rod_frac)     # set lever position from inserted fraction
        dial_detents[0] = ch2_det                     # set CH2 dial detent
        dial_detents[1] = ch3_det                     # set CH3 dial detent

        if power_sw and not power_on:                 # rising edge: switch OFF → ON from Arduino telemetry
            power_on = True                           # enable power state
            scram_active = False                      # clear scram state
            scram_start_time = 0.0                    # reset scram timing
            scram_reason = ""                         # clear scram reason
            reset_on_power_on(now)                    # reset physics and instruments to startup conditions
        if (not power_sw) and power_on and (not scram_active):  # falling edge: ON → OFF triggers scram
            start_scram("POWER SWITCH OFF", now)       # enter SCRAM sequence
        if scram_sw and power_on and (not scram_active):         # external scram input
            start_scram("EXTERNAL SCRAM", now)         # enter SCRAM sequence

    # ---------- EVENTS ----------
    for event in pygame.event.get():                  # process Pygame event queue (for windowed mode)
        if event.type == pygame.QUIT:                 # window close button (for windowed mode)
            running = False                           # exit loop (for windowed mode)

        elif event.type == pygame.KEYDOWN:            # ESC key exits fullscreen
            if event.key == pygame.K_ESCAPE:          # window close execute
                running = False                       # exit loop

        elif event.type == pygame.MOUSEBUTTONDOWN:    # mouse click
            mx, my = event.pos                        # mouse coordinates (px)

            if TOGGLE_RECT.collidepoint((mx, my)):    # if click on “Show/Hide Controls” tab
                controls_collapsed = not controls_collapsed  # toggle state
                update_layout(controls_collapsed)     # recompute layout
                selected_lever = False                # drop selections when collapsing
                lever_active = False
                selected_dial = None
                dial_active[:] = [False, False]
                continue                              # done handling click

            if controls_collapsed:                    # if controls hidden, ignore clicks on control region
                continue

            if POWER_RECT.collidepoint((mx, my)):     # if click POWER button
                if transport is not None:             # when connected, send the exact toggle command expected by the Arduino sketch
                    send_flipmask_command(transport, 1 << 1)
                else:                                 # fallback standalone behavior when no Arduino is attached
                    if power_on:                      # if currently on, turn off via scram sequence
                        start_scram("POWER BUTTON OFF", now)
                    else:                             # if currently off, start the system
                        power_on = True
                        scram_active = False
                        scram_start_time = 0.0
                        scram_reason = ""
                        reset_on_power_on(now)
                continue

            if SCRAM_RECT.collidepoint((mx, my)):     # if click SCRAM button
                if power_on and not scram_active:     # only scram if running and not already scramming
                    if transport is not None:         # send Arduino command when hardware is attached
                        send_flipmask_command(transport, 1 << 0)
                    else:                             # fallback standalone behavior with no Arduino connected
                        start_scram("MANUAL SCRAM", now)
                continue

            if power_on and not scram_active:         # only allow selecting controls when running normally
                if (LEVER_X - 14 < mx < LEVER_X + 14 and lever_y < my < lever_y + LEVER_HEIGHT):  # lever hitbox
                    selected_lever = True             # mark lever selected
                    lever_active = True               # highlight lever active
                    selected_dial = None              # deselect any dial
                    dial_active[:] = [False, False]   # deactivate both dials
                    lever_activation_time = now       # update activity timestamp
                else:
                    selected_lever = False            # clicking elsewhere deselects lever
                    lever_active = False

                for i, (cx, cy) in enumerate(dial_centers):     # check dial hit tests
                    if math.hypot(mx - cx, my - cy) < DIAL_RADIUS + 12:  # within dial radius
                        selected_dial = i             # select that dial index
                        dial_active[:] = [False, False]
                        dial_active[i] = True         # activate selected dial highlight
                        selected_lever = False        # deselect lever
                        lever_active = False
                        dial_activation_time[i] = now # update dial activity timestamp
                        break                         # stop after selecting one dial

    # ---------- KEYBOARD CONTROL ----------
    keys = pygame.key.get_pressed()                   # current keyboard state snapshot
    if (not controls_collapsed) and power_on and not scram_active:  # only allow keyboard control when visible & running
        if selected_lever and lever_active:           # if lever is selected and active
            speed = 260.0                             # lever speed (px/s) when holding arrow key
            if keys[pygame.K_UP]:                     # up arrow inserts rod (move handle up)
                lever_y -= speed * dt                 # update lever y by speed*dt
                lever_activation_time = now           # refresh activity timer
            elif keys[pygame.K_DOWN]:                 # down arrow withdraws rod (move handle down)
                lever_y += speed * dt
                lever_activation_time = now
            lever_y = clamp(lever_y, LEVER_TRACK_TOP, LEVER_TRACK_BOTTOM - LEVER_HEIGHT)  # keep within track

        if selected_dial is not None and dial_active[selected_dial]: # if a dial is selected and active
            if keys[pygame.K_LEFT]:                   # left arrow decreases detent
                dial_detents[selected_dial] = clamp(dial_detents[selected_dial] - 1, -3, 3)
                dial_activation_time[selected_dial] = now
            elif keys[pygame.K_RIGHT]:                # right arrow increases detent
                dial_detents[selected_dial] = clamp(dial_detents[selected_dial] + 1, -3, 3)
                dial_activation_time[selected_dial] = now

    # ---------- COOLDOWN AUTO-DESELECT ----------
    if lever_active and (now - lever_activation_time) > COOLDOWN_TIME:  # if lever idle too long
        lever_active = False                     # deactivate lever highlight
        selected_lever = False                   # deselect lever
    for i in range(NUM_DIALS):                   # for each dial
        if dial_active[i] and (now - dial_activation_time[i]) > COOLDOWN_TIME:  # if idle too long
            dial_active[i] = False               # deactivate dial highlight
            if selected_dial == i:               # if that dial was selected
                selected_dial = None             # clear selection

    # ---------- SCRAM BEHAVIOR ----------
    if scram_active:                             # if currently scramming
        lever_y = LEVER_TRACK_TOP                # force lever fully inserted visually
        if (now - scram_start_time) > 3.0:       # after scram hold duration
            power_on = False                     # shut system power off
            scram_active = False                 # exit scram state
            scram_start_time = 0.0               # clear scram timing
            scram_reason = ""                    # clear reason

            power_history.clear()                # clear plot histories
            ch2_history.clear()
            ch3_history.clear()

            n = 1e-12                            # reset neutron state to floor
            C = 0.0                              # reset precursor
            rho_prev = 0.0                       # reset reactivity
            rho_now  = 0.0
            ch2_I_ind = I_FLOOR                  # reset instrument indications
            ch3_I_ind = I_FLOOR

            sim_time_s = 0.0                     # reset sim time
            power_on_wall_time = 0.0

            high_scram_timer = 0.0               # reset debounce timers
            low_scram_timer  = 0.0

            lever_y = LEVER_TRACK_TOP            # keep lever inserted
            dial_detents[:] = [0, 0]             # reset dials to center detent
            lever_active = False                 # drop selections
            dial_active[:] = [False, False]
            selected_dial = None
            selected_lever = False

    # ---------- PHYSICS ----------
    power_raw = 0.0                              # default power if not running
    if power_on and not scram_active:            # only update physics when system is on and not in scram
        sim_time_s += dt                         # advance simulation time

        x_inserted = lever_inserted_fraction(lever_y)        # convert lever y → inserted fraction 0..1
        f_inserted = fcr_inserted_worth_fraction(x_inserted) # curve mapping: fraction of worth inserted

        rho_cmd = RHO_WITHDRAWN + RHO_WORTH_MAG * f_inserted # commanded reactivity from rod worth curve
        rho_cmd = clamp(rho_cmd, RHO_NEG_LIMIT, RHO_POS_LIMIT) # enforce safety bounds

        if abs(rho_cmd - rho_prev) < RHO_DEADBAND:           # deadband to prevent chatter
            rho_cmd = rho_prev

        rho_old = rho_prev                                   # store previous applied reactivity for prompt-jump logic
        rho_now = first_order_lag(rho_prev, rho_cmd, dt, RHO_TAU)  # apply lag filter to reactivity

        if PROMPT_JUMP_ENABLE and abs(rho_now - rho_old) > PROMPT_JUMP_EPS:  # apply prompt jump only for meaningful Δρ
            n = apply_prompt_jump(n, rho_old, rho_now)       # instantaneous neutron-population jump; C does not jump

        substeps = 12                           # subdivide dt to keep RK4 stable for small Λ
        h = dt / substeps                       # substep size
        for _ in range(substeps):              # integrate across dt using many small steps
            n, C = rk4_step(n, C, h, rho_now)  # advance kinetics states after the prompt jump
            if (not (n == n)) or n < 1e-12:    # guard NaN and negative state
                n = 1e-12
            if (not (C == C)) or C < 0.0:      # guard NaN and negative precursor
                C = 0.0

        rho_prev = rho_now                     # store applied ρ for next frame
        power_raw = max(0.0, n)               # displayed power is taken directly from kinetics amplitude

        # ---------- AUTO-SCRAM LOGIC ----------
        if power_raw >= SCRAM_HIGH_W:          # if above high threshold
            high_scram_timer += dt             # accumulate time above threshold
        else:
            high_scram_timer = 0.0             # reset if condition clears

        if power_raw <= SCRAM_LOW_W:           # if below low threshold
            low_scram_timer += dt              # accumulate time below threshold
        else:
            low_scram_timer = 0.0              # reset if condition clears

        if high_scram_timer >= SCRAM_DEBOUNCE_S:  # if high threshold sustained long enough
            start_scram("AUTO SCRAM: HIGH POWER", now)

        if low_scram_timer >= SCRAM_DEBOUNCE_S:   # if low threshold sustained long enough
            start_scram("AUTO SCRAM: ZERO POWER", now)

        power_plot = clamp(power_raw, 0.0, MAX_POWER_DISPLAY)  # cap only for plot display
        power_history.append(power_plot)       # store for history buffer

        ch2_I_true = max(I_FLOOR, I_PER_WATT_CH2 * power_raw)  # “true” current for CH2 based on power
        ch3_I_true = max(I_FLOOR, I_PER_WATT_CH3 * power_raw)  # “true” current for CH3 based on power

        ch2_I_ind = first_order_lag(ch2_I_ind, ch2_I_true, dt, CH2_TAU)  # instrument lag on CH2
        ch3_I_ind = first_order_lag(ch3_I_ind, ch3_I_true, dt, CH3_TAU)  # instrument lag on CH3

        ch2_history.append(ch2_I_ind)          # push indicated current sample into CH2 trace
        ch3_history.append(ch3_I_ind)          # push indicated current sample into CH3 trace
    else:
        ch2_I_ind = first_order_lag(ch2_I_ind, I_FLOOR, dt, 0.5)  # decay indication toward floor when off/scram
        ch3_I_ind = first_order_lag(ch3_I_ind, I_FLOOR, dt, 0.5)

    # ---------- BLUETOOTH DATA STREAM ----------
    # Sends current simulation time and neutron state and power to the second Pi without ever stalling the GUI
    # Try to connect only occasionally (NOT every frame), and never block the loop for long
    if BT_ENABLE and bt_sock is None:
        if now >= bt_enable_time and now >= bt_next_retry: # only attempt connection at a controlled rate
            bt_next_retry = now + BT_RETRY_S              # schedule the next allowed retry time
            try:
                bt_sock = bt_connect()                    # attempt to connect (has short timeout)
            except Exception:
                bt_sock = None                            # stay disconnected; try again in soon

    if BT_ENABLE and bt_sock is not None:
        if (now - bt_last_send) >= (1.0 / BT_HZ):         # If connected, send at BT_HZ and send at a fixed rate
            bt_last_send = now
            try:
                line = f"{sim_time_s:.4f},{n:.6e}\n"      # Format: sim_time_s, n
                bt_sock.send(line.encode("utf-8"))        # transmit
            except (BlockingIOError, InterruptedError):
                pass                                      # skip this frame if non-blocking Bluetooth send would stall
            except Exception:
                try:
                    bt_sock.close()                       # close broken socket
                except Exception:
                    pass
                bt_sock = None                            # force reconnect path next time
                bt_next_retry = now + BT_RETRY_S          # back off before retrying

    # =========================
    # DRAW PLOTS
    # =========================
    pygame.draw.rect(screen, GREY, (GRAPH_X, GRAPH_Y, GRAPH_WIDTH, GRAPH_HEIGHT))          # plot background
    pygame.draw.rect(screen, BLACK, (GRAPH_X, GRAPH_Y, linear_width, GRAPH_HEIGHT), 2)    # left plot border
    pygame.draw.rect(screen, BLACK, (right_start, GRAPH_Y, linear_width, GRAPH_HEIGHT), 2)# right plot border
    pygame.draw.line(screen, BLACK, (right_start, GRAPH_Y), (right_start, GRAPH_Y + GRAPH_HEIGHT), 2)  # divider line

    # CH2 LOG (left)
    if len(ch2_history) > 1:                   # only draw if there are enough points
        det2 = dial_detents[0]                 # CH2 dial detent value
        ch2_log_max = CH2_LOG_MAX_AT_0 + det2  # top decade shifts with detent
        ch2_log_min = ch2_log_max - CH2_LOG_SPAN_DECADES  # bottom decade is fixed span below top
        pts = []                               # list of polyline points
        L = len(ch2_history)                   # number of samples
        for i, I in enumerate(ch2_history):    # iterate history
            y = GRAPH_Y + (1.0 - i / max(1, L - 1)) * GRAPH_HEIGHT  # map sample index → vertical position (time axis)
            logI = math.log10(max(I, I_FLOOR)) # convert current to log10, clamp to avoid log(0)
            frac = (logI - ch2_log_min) / max(1e-9, (ch2_log_max - ch2_log_min)) # normalize into [0,1] window
            frac = clamp(frac, 0.0, 1.0)       # clamp inside the window
            x = GRAPH_X + frac * linear_width  # map normalized fraction → horizontal position
            pts.append((x, y))                 # append point
        pygame.draw.lines(screen, YELLOW, False, pts, 2)  # draw log trace

    # CH3 LIN (right)
    if len(ch3_history) > 1:                   # only draw if there are have enough points
        det3 = dial_detents[1]                 # CH3 dial detent value
        fs3 = CH3_RANGES[det3 + 3]             # map detent -3..3 → list index 0..6
        pts = []                               # list of polyline points
        L = len(ch3_history)                   # number of samples
        for i, I in enumerate(ch3_history):    # iterate history
            y = GRAPH_Y + (1.0 - i / max(1, L - 1)) * GRAPH_HEIGHT  # map sample index → vertical position
            frac = clamp(I / max(fs3, I_FLOOR), 0.0, 1.0)     # linear fraction of full-scale
            x = right_start + frac * linear_width             # map fraction → horizontal position
            pts.append((x, y))                 # append point
        pygame.draw.lines(screen, RED, False, smooth_points(pts, steps=4), 2)  # draw smoothed linear trace with lighter densification for less per-frame work

    # =========================
    # HEADER TEXT
    # =========================
    header_y0 = GRAPH_TOP + 10                 # header baseline y position

    n_line = inst_font.render(                 # render neutron state and sim time
        f"Neutron state n: {n:.3e}   Sim time: {sim_time_s:8.2f} s",
        True, BLACK
    )
    screen.blit(n_line, (GRAPH_X, header_y0))  # draw the neutron/time line left-aligned

    transport_age = (now - transport_last_rx_time) if transport_last_rx_time > 0.0 else 1e9  # age of latest good Arduino packet
    transport_state = "LIVE" if transport_age < 1.0 else ("CONNECTED" if transport is not None else "SEARCHING")
    transport_line = inst_font.render(
        f"Link: {transport_name}   mode: {TRANSPORT_MODE}   state: {transport_state}   RX: {transport_rx_packets}   TX: {transport_tx_commands}",
        True, BLACK
    )
    screen.blit(transport_line, (GRAPH_X, header_y0 + 18))

    k_eff = 1.0 / max(1e-9, (1.0 - rho_now))   # approximate k_eff from ρ using ρ=(k-1)/k => k=1/(1-ρ)
    dollars = rho_now / beta                   # reactivity in dollars ($ = ρ/β)

    p_text = font.render(                      # render power readout
        f"Power (raw): {format_watts(power_raw)}",
        True, BLACK
    )
    screen.blit(p_text, p_text.get_rect(center=(GRAPH_X + GRAPH_WIDTH // 2, header_y0 + 12)))  # centered

    k_text = inst_font.render(                 # render k_eff, rho, dollars line
        f"k_eff: {k_eff:.6f}   ρ: {rho_now:+.6f}   $: {dollars:+.2f}",
        True, BLACK
    )
    screen.blit(k_text, k_text.get_rect(center=(GRAPH_X + GRAPH_WIDTH // 2, header_y0 + 52)))  # centered

    if transport_last_error:
        err = inst_font.render(f"Link error: {transport_last_error}", True, BLACK)
        screen.blit(err, (GRAPH_X, header_y0 + 36))

    left_label = title_font.render("Channel 2 (Log Current)", True, BLACK)  # left plot title
    screen.blit(left_label, left_label.get_rect(center=(GRAPH_X + linear_width // 2, GRAPH_Y - 18)))  # place above left plot
    right_label = title_font.render("Channel 3 (Linear Current)", True, BLACK) # right plot title
    screen.blit(right_label, right_label.get_rect(center=(right_start + linear_width // 2, GRAPH_Y - 18))) # place above right plot

    if scram_reason:                            # if scram reason exists, display it under graphs
        sr = inst_font.render(f"SCRAM Reason: {scram_reason}", True, BLACK)
        screen.blit(sr, (GRAPH_X, GRAPH_Y + GRAPH_HEIGHT + TOGGLE_GAP + TOGGLE_H + 6))

    # =========================
    # COLLAPSE/EXPAND TAB (always visible)
    # =========================
    pygame.draw.rect(screen, (220, 220, 220), TOGGLE_RECT)  # draw tab background
    pygame.draw.rect(screen, BLACK, TOGGLE_RECT, 2)         # draw tab border
    tab_txt = "▼ Show Controls" if controls_collapsed else "▲ Hide Controls"  # tab label depends on state
    tab = tab_font.render(tab_txt, True, BLACK)             # render tab label
    screen.blit(tab, tab.get_rect(center=TOGGLE_RECT.center)) # center label in tab

    # =========================
    # CONTROLS (only if expanded)
    # =========================
    if not controls_collapsed:                   # only draw and interact with controls if visible
        pygame.draw.rect(screen, GREEN if power_on else GREY, POWER_RECT)  # draw POWER button
        btn = font.render("POWER", True, BLACK)  # render “POWER” text
        screen.blit(btn, btn.get_rect(center=POWER_RECT.center))  # center on button

        pygame.draw.rect(screen, SCRAM_RED if scram_active else GREY, SCRAM_RECT) # draw SCRAM button
        scr = font.render("SCRAM", True, BLACK)   # render “SCRAM” text
        screen.blit(scr, scr.get_rect(center=SCRAM_RECT.center)) # center label inside SCRAM button
        pygame.draw.rect(screen, GREY, (LEVER_X - 6, LEVER_TRACK_TOP, 12, LEVER_TRACK_BOTTOM - LEVER_TRACK_TOP)) # lever track
        handle_color = LEVER_INACTIVE_COLOR if (scram_active or not power_on) else (YELLOW if lever_active else BLUE) # handle state color
        pygame.draw.rect(screen, handle_color, (LEVER_X - 14, lever_y, 28, LEVER_HEIGHT)) # lever handle
        nm = font.render(lever_name, True, BLACK) # lever name label
        screen.blit(nm, nm.get_rect(center=(LEVER_X, LEVER_TRACK_BOTTOM + 34))) # place under lever

        for i in range(NUM_DIALS):                 # draw each dial
            cx, cy = dial_centers[i]               # dial center
            pygame.draw.circle(screen, GREY, (cx, cy), DIAL_RADIUS + 8) # outer ring
            face = DIAL_INACTIVE_COLOR if (not power_on) else (YELLOW if dial_active[i] else DIAL_INACTIVE_COLOR) # face color
            pygame.draw.circle(screen, face, (cx, cy), DIAL_RADIUS)     # dial face
            draw_dial_ticks(cx, cy, DIAL_RADIUS)   # tick marks and labels

            ang = math.radians(angle_from_detent(dial_detents[i])) # detent → degrees → radians
            x2 = cx + DIAL_RADIUS * math.cos(ang)  # pointer endpoint x
            y2 = cy + DIAL_RADIUS * math.sin(ang)  # pointer endpoint y
            pygame.draw.line(screen, WHITE, (cx, cy), (x2, y2), 4) # pointer line

            cap = font.render(dial_names[i], True, BLACK)          # dial name label
            screen.blit(cap, cap.get_rect(center=(cx, cy - DIAL_RADIUS - 34))) # place above dial

            det_txt = inst_font.render(f"det {dial_detents[i]:+d}", True, BLACK) # detent readout
            screen.blit(det_txt, det_txt.get_rect(center=(cx, cy + DIAL_RADIUS + 14))) # place below dial

    pygame.display.flip()           # present the composed frame to the screen
    clock.tick(60)                  # cap loop to ~60 frames per second (also stabilizes dt)

pygame.quit()                      # shut down Pygame cleanly
sys.exit()                         # exit program
