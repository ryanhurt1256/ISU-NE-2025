import pygame        # pygame: window, drawing, events, timing
import sys           # sys: used for sys.exit at the end
import math          # math: exp, log10, trig functions, etc.
import time          # time: used to compute frame-to-frame time step dt
import collections   # collections: deque for rolling power history

try:
    import serial    # try to import pyserial for ESP32 communications
except Exception:
    serial = None    # if import fails, set serial to None to disable serial I/O

pygame.init()        # initialize all pygame subsystems (video, font, etc.)

WIDTH, HEIGHT = 1000, 1000                                         # window width and height in pixels
screen = pygame.display.set_mode((WIDTH, HEIGHT))                  # main drawing surface
pygame.display.set_caption("Control Panel")                        # title shown on the window bar
clock = pygame.time.Clock()                                        # pygame clock to cap FPS and get dt

WHITE  = (255, 255, 255)                                           # color constant: white (background)
GREY   = (150, 150, 150)                                           # color constant: grey (graph fill, lever track)
BLUE   = (50, 100, 255)                                            # color constant: blue (inactive lever handle)
RED    = (255, 60, 60)                                             # color constant: red (linear power curve)
YELLOW = (255, 255, 0)                                             # color constant: yellow (log power curve)
BLACK  = (0, 0, 0)                                                 # color constant: black (axes, text)
SCRAM_RED = (255, 0, 0)                                            # color constant: bright red (SCRAM button)
DIAL_INACTIVE_COLOR  = (50, 50, 50)                                # color constant: dark grey dial when idle
LEVER_INACTIVE_COLOR = (80, 80, 80)                                # color constant: lever when power off / scram
GREEN  = (0, 200, 0)                                               # color constant: green (POWER on button)

MARGIN_LEFT, MARGIN_TOP, MARGIN_RIGHT, MARGIN_BOTTOM = 140, 70, 100, 150  # outer margins around the graph
GRAPH_WIDTH  = WIDTH - MARGIN_LEFT - MARGIN_RIGHT                 # width available for both linear+log panels
GRAPH_HEIGHT = 300                                                # fixed vertical size of graph in pixels
GRAPH_X, GRAPH_Y = MARGIN_LEFT, MARGIN_TOP                        # top-left corner of graph area

LEVER_TRACK_TOP    = GRAPH_Y + GRAPH_HEIGHT + 110                 # y-position at which lever handles start (fully out)
LEVER_TRACK_BOTTOM = LEVER_TRACK_TOP + 220                        # y-position for bottom of lever travel (fully in)
LEVER_HEIGHT       = 80                                           # height of the drawn lever handle rectangle
NUM_LEVERS         = 4                                            # number of control rods / levers
lever_spacing      = GRAPH_WIDTH // (NUM_LEVERS + 1)              # equal horizontal spacing between lever centers
lever_x_positions  = [GRAPH_X + lever_spacing * (i + 1)           # compute x-position for each lever center
                      for i in range(NUM_LEVERS)]
lever_positions    = [LEVER_TRACK_TOP] * NUM_LEVERS               # all rods start visually fully out (top position)

NUM_DIALS   = 2                                                   # two dials: one for linear, one for log zoom (visual)
DIAL_RADIUS = 40                                                  # radius of drawn dials in pixels
dial_spacing = GRAPH_WIDTH // (NUM_DIALS + 1)                     # horizontal spacing for dial centers
dial_centers = [(GRAPH_X + dial_spacing * (i + 1),                # centers of dials at same y below levers
                 LEVER_TRACK_BOTTOM + 200)
                for i in range(NUM_DIALS)]
dial_angles  = [-90] * NUM_DIALS                                  # each dial starts at -90°, which we treat as "0" mark

selected_lever = None                                             # index of lever currently selected with mouse, else None
selected_dial  = None                                             # index of dial currently selected with mouse, else None
lever_active   = [False] * NUM_LEVERS                             # True while a lever is active for keyboard movement
dial_active    = [False] * NUM_DIALS                              # True while a dial is active for keyboard rotation
lever_activation_time = [0] * NUM_LEVERS                          # last time stamp of activity per lever
dial_activation_time  = [0] * NUM_DIALS                           # last time stamp of activity per dial
COOLDOWN_TIME = 3                                                 # if no activity for 3 seconds → auto-deselect control

scram_active     = False                                          # True while SCRAM state is in effect
scram_start_time = 0                                              # time at which SCRAM was initiated

power_history = collections.deque(maxlen=300)                     # fixed-length buffer of power samples (newest at end)

MAX_POWER = 5.0                                                   # linear panel range: 0–5 W full-scale
P_SCALE   = 3.0                                                   # mapping from neutron population n to displayed power P
current_power = 3.0                                               # initial power: 3 W = mid-range on linear panel

LOG_MIN   = -3.0                                                  # left log-panel edge: 10^(-3) W
LOG_MAX   =  1.0                                                  # right log-panel edge: 10^(1) W = 10 W
P_MIN_LOG = 10 ** LOG_MIN                                         # minimum power used when taking log10 (avoid log(0))

beta   = 0.0065                                                   # effective delayed neutron fraction β (dimensionless)
Lambda = 5e-5                                                     # prompt neutron lifetime Λ in seconds (very small)
lam    = 0.08                                                     # effective single delayed neutron group decay constant λ

S_SOURCE = 201                                                    # external neutron source term (normalized units)

n = current_power / P_SCALE                                       # neutron population n inferred from starting power P
C = (beta / (Lambda * lam)) * n                                   # precursor concentration C at steady state: dC/dt=0
rho_prev = 0.0                                                    # previous reactivity value used for prompt jump calc

last_time = time.time()                                           # last frame time stamp for dt computation

button_width, button_height = 80, 40                              # POWER button dimensions
power_on = False                                                  # whether reactor is considered "running"
power_button_rect = pygame.Rect(60, 875, button_width, button_height)  # rectangle defining POWER button hitbox

lever_names = ["Fine CR", "Fine CR", "Fine CR", "Coarse CR"]   # labels printed below each lever
dial_names  = ["Lin Zoom", "Log Zoom"]                            # labels printed above each dial

def rod_curve(x: float) -> float:
    """Return fractional worth f(x) of a rod at fractional insertion x."""
    a = 3.3                                                       # exponential falloff constant fitted to AGN data
    return 1.0 - math.exp(-a * x)                                # f(x) = 1 − e^(−ax), 0 at x=0, approx 1 at x=1

def smooth_points(points, steps=5):
    """Insert intermediate points between each pair to make smoother plotted curves."""
    smoothed = []                                                 # list to accumulate new interpolated points
    if len(points) < 2:                                          # if fewer than 2 points, nothing to smooth
        return points
    for i in range(len(points) - 1):                             # iterate over segments between point i and i+1
        x1, y1 = points[i]                                       # starting point of segment
        x2, y2 = points[i + 1]                                   # ending point of segment
        for t in range(steps):                                   # create 'steps' intermediate points
            alpha = t / steps                                    # normalized interpolation factor in [0,1)
            sx = x1 * (1 - alpha) + x2 * alpha                   # linear interpolation for x-coordinate
            sy = y1 * (1 - alpha) + y2 * alpha                   # linear interpolation for y-coordinate
            smoothed.append((sx, sy))                            # store interpolated point
    smoothed.append(points[-1])                                  # include final original endpoint
    return smoothed                                              # return densified polyline

def draw_ticks_and_labels():
    """Draw bottom x-axis (linear+log) and left y-axis tick marks + labels."""
    font = pygame.font.SysFont(None, 20)                         # font for axis labels, modest size
    linear_width = GRAPH_WIDTH // 2                              # width occupied by linear panel
    right_start  = GRAPH_X + linear_width                        # x location where log panel starts

    for power_tick in range(0, int(MAX_POWER) + 1):              # loop over integer powers 0..5 W
        x = GRAPH_X + (power_tick / MAX_POWER) * linear_width    # map tick position into [GRAPH_X, GRAPH_X+linear_width]
        y = GRAPH_Y + GRAPH_HEIGHT                               # y location at bottom of graph
        pygame.draw.line(screen, BLACK, (x, y), (x, y + 8), 2)   # draw short vertical tick mark
        if power_tick == int(MAX_POWER):                         # skip label at 5 W to avoid overlap at split
            continue
        label = font.render(str(power_tick), True, BLACK)        # render tick text (e.g., "0","1",...)
        screen.blit(label, label.get_rect(center=(x, y + 20)))    # center label under tick

    log_ticks  = [LOG_MIN, -2.0, -1.0, 0.0, 1.0]                 # log10 positions for 10^-3,10^-2,10^-1,10^0,10^1
    log_labels = [" ", "10^-2", "10^-1", "10^0", "10^1"]         # first blank so there is no label at split
    for tick, lbl in zip(log_ticks, log_labels):                 # pair each tick with its label
        frac = (tick - LOG_MIN) / (LOG_MAX - LOG_MIN)            # normalize tick position to [0,1] across log range
        frac = max(0.0, min(1.0, frac))                          # clamp to [0,1] to avoid rounding overshoot
        x = right_start + frac * linear_width                    # map fraction into log panel x range
        y = GRAPH_Y + GRAPH_HEIGHT                               # y at bottom of graph
        pygame.draw.line(screen, BLACK, (x, y), (x, y + 8), 2)   # draw log tick
        if lbl.strip() == "":                                    # if label is blank string, do not draw any text
            continue
        text = font.render(lbl, True, BLACK)                     # render log tick label (10^-2, etc.)
        screen.blit(text, text.get_rect(center=(x, y + 20)))     # place label below tick

    max_time, tick_step = power_history.maxlen, 50               # max_time is length of deque (300); step 50 for ticks
    for time_tick in range(0, max_time + 1, tick_step):          # ticks at 0,50,100,...,300
        y = GRAPH_Y + (time_tick / max_time) * GRAPH_HEIGHT      # map time index to vertical position in graph
        x = GRAPH_X                                              # left edge of graph
        pygame.draw.line(screen, BLACK, (x - 8, y), (x, y), 2)   # small horizontal tick on left margin
        label = font.render(str(time_tick), True, BLACK)         # label with index number (proxy for time)
        screen.blit(label, label.get_rect(center=(x - 25, y)))   # position labels left of ticks

def draw_dial_ticks(cx, cy, radius):
    """Draw static tick marks and -3..+3 labels around a dial centered at (cx,cy)."""
    font = pygame.font.SysFont(None, 18)                         # slightly smaller font than axes
    num_marks, angle_step = 7, 30                                # 7 marks, spaced 30° around dial
    for i in range(num_marks):                                   # i = 0..6, center mark is i=3
        angle_deg = -90 - angle_step * (i - 3)                   # main 'zero' tick at -90°, others symmetric
        angle_rad = math.radians(angle_deg)                      # convert degrees to radians for trig
        tick_len  = 10 if i == 3 else 5                          # center tick longer so it stands out
        label_txt = "0" if i == 3 else str(i - 3)                # center label "0", others "-3,-2,-1,1,2,3"
        outer_x = cx + radius * math.cos(angle_rad)              # x of outer end of tick on circle
        outer_y = cy + radius * math.sin(angle_rad)              # y of outer end of tick on circle
        inner_x = cx + (radius - tick_len) * math.cos(angle_rad) # x of inner end (slightly inwards)
        inner_y = cy + (radius - tick_len) * math.sin(angle_rad) # y of inner end
        pygame.draw.line(screen, BLACK, (outer_x, outer_y), (inner_x, inner_y), 2)  # draw tick line
        label_x = cx + (radius + 22) * math.cos(angle_rad)       # x for label outside tick
        label_y = cy + (radius + 22) * math.sin(angle_rad)       # y for label
        label = font.render(label_txt, True, BLACK)              # render dial label text
        screen.blit(label, label.get_rect(center=(label_x, label_y)) )  # draw dial label

SERIAL_PORT, BAUDRATE = "/dev/ttyUSB0", 115200                   # default serial port and baud rate for ESP32

def setup_serial(port=SERIAL_PORT, baud=BAUDRATE):
    """Try to open a serial port and return the Serial object or None on failure."""
    if serial is None:                                           # pyserial is not available
        print("pyserial not installed or unavailable; serial input disabled.")
        return None
    try:
        ser = serial.Serial(port, baud, timeout=0.05)            # open serial with 50 ms timeout
        time.sleep(1.0)                                          # give MCU a moment to reset
        print("Serial connected on", port)                       # debug message
        return ser
    except Exception as e:                                       # any error while opening serial
        print("Serial connect failed:", e)                       # print error message
        return None                                              # disable serial support

ser = setup_serial()                                             # create serial connection once at startup

def read_controls_from_serial(ser):
    """Read one line from serial and parse rod fractions + switches; return tuple or None."""
    if ser is None:                                              # if serial is disabled, nothing to read
        return None
    try:
        raw = ser.readline().decode(errors="ignore").strip()     # read text line, remove trailing whitespace
        if not raw:                                              # if empty string → timeout or nothing received
            return None
        parts = raw.split(",")                                   # split comma-separated values
        if len(parts) < 6:                                       # need at least 6 values: 4 rods + power + scram
            return None
        vals     = [float(p) for p in parts[:4]]                 # rods represented as floats in [0,1]
        power_sw = bool(int(float(parts[4])))                    # convert 5th element to bool: power switch on/off
        scram_sw = bool(int(float(parts[5])))                    # convert 6th to bool: scram switch on/off
        return vals, power_sw, scram_sw                          # return parsed control tuple
    except Exception:                                            # swallow any parse/serial error
        return None

def kinetics_rhs(state, rho):
    """Compute derivatives (dn/dt, dC/dt) for 1-group point kinetics with source."""
    n_s, C_s = state                                             # unpack current neutron population and precursors
    dn_dt = ((rho - beta) / Lambda) * n_s + lam * C_s + S_SOURCE # dn/dt: prompt term + delayed term + external source
    dC_dt = (beta / Lambda) * n_s - lam * C_s                    # dC/dt: production from fission minus decay
    return dn_dt, dC_dt                                          # return time derivatives

def rk4_step(state, dt_sub, rho):
    """Advance state (n,C) over time dt_sub using classical 4th-order Runge–Kutta for ODEs."""
    n0, C0 = state                                               # starting values at t
    k1n, k1C = kinetics_rhs((n0, C0), rho)                       # slope at beginning of interval
    k2n, k2C = kinetics_rhs((n0 + 0.5 * dt_sub * k1n,           # slope at midpoint using k1
                             C0 + 0.5 * dt_sub * k1C), rho)
    k3n, k3C = kinetics_rhs((n0 + 0.5 * dt_sub * k2n,           # slope at midpoint using k2
                             C0 + 0.5 * dt_sub * k2C), rho)
    k4n, k4C = kinetics_rhs((n0 + dt_sub * k3n,                 # slope at end using k3
                             C0 + dt_sub * k3C), rho)
    n_new = n0 + (dt_sub / 6.0) * (k1n + 2 * k2n + 2 * k3n + k4n)  # RK4 weighted average for n
    C_new = C0 + (dt_sub / 6.0) * (k1C + 2 * k2C + 2 * k3C + k4C)  # RK4 weighted average for C
    return n_new, C_new                                          # return next-step (n,C)

running = True                                                   # main loop flag
while running:                                                   # loop until user quits
    screen.fill(WHITE)                                           # clear screen each frame
    now = time.time()                                            # current wall clock time in seconds
    dt = now - last_time                                         # raw frame time step since last frame
    if dt <= 0:                                                  # if dt is zero or negative (rare clock issue)
        dt = 1.0 / 60.0                                          # fall back to ~16.7 ms (60 FPS)
    if dt > 0.1:                                                 # if dt is very large (window dragged, etc.)
        dt = 0.1                                                 # clamp dt to avoid huge physics jumps
    last_time = now                                              # update last_time for next frame

    serial_controls = read_controls_from_serial(ser)             # attempt to read ESP32 control line
    if serial_controls:                                          # if valid data
        rods, power_sw, scram_sw = serial_controls               # unpack rod fractions and switches
        travel = (LEVER_TRACK_BOTTOM - LEVER_TRACK_TOP - LEVER_HEIGHT)  # maximum vertical movement for handles
        for i in range(min(NUM_LEVERS, len(rods))):              # update lever positions for each rod
            frac = max(0.0, min(1.0, rods[i]))                   # clamp rod fraction into [0,1]
            lever_positions[i] = LEVER_TRACK_TOP + frac * travel # map fractional insertion → screen coordinate
        if power_sw and not power_on:                            # hardware power switch turned ON
            power_on = True                                      # mark reactor as running
            current_power = 3.0                                  # reset indicated power to 3 W
            n = current_power / P_SCALE                          # recompute neutron population from power
            C = (beta / (Lambda * lam)) * n                      # set precursors to steady-state level for this n
            power_history.clear()                                # clear power history graph
        elif not power_sw and power_on:                          # hardware power switch turned OFF
            scram_active = True                                  # treat as SCRAM
            scram_start_time = now                               # remember SCRAM start time
        if scram_sw:                                             # dedicated SCRAM switch pressed
            scram_active = True                                  # enter SCRAM
            scram_start_time = now                               # record time

    if scram_active and (now - scram_start_time > 3):            # if SCRAM has lasted more than 3 seconds
        power_on = False                                         # turn power state off
        scram_active = False                                     # end SCRAM state
        current_power = 3.0                                      # reset indicated power back to 3 W
        n = current_power / P_SCALE                              # recompute neutron population from power
        C = (beta / (Lambda * lam)) * n                          # reset precursors to steady-state
        lever_positions[:] = [LEVER_TRACK_TOP] * NUM_LEVERS      # move all levers back to fully out
        dial_angles[:]     = [-90] * NUM_DIALS                   # reset dial angles to neutral
        selected_lever = selected_dial = None                    # clear active selections
        lever_active[:] = [False] * NUM_LEVERS                   # mark all levers inactive
        dial_active[:]  = [False] * NUM_DIALS                    # mark all dials inactive
        power_history.clear()                                    # clear power history

    for event in pygame.event.get():                             # process all pending pygame events
        if event.type == pygame.QUIT:                            # window close button pressed
            running = False                                      # break main loop on next iteration
        elif event.type == pygame.MOUSEBUTTONDOWN:               # mouse button pressed
            mx, my = event.pos                                   # mouse coordinates
            if power_button_rect.collidepoint((mx, my)):         # click occurred inside POWER button
                if power_on:                                     # if power currently on
                    scram_active = True                          # initiate SCRAM sequence
                    scram_start_time = now                       # record start time
                else:                                            # power currently off
                    power_on = True                              # turn reactor on
                    current_power = 3.0                          # set initial power 3 W
                    n = current_power / P_SCALE                  # update neutron population
                    C = (beta / (Lambda * lam)) * n              # update precursors
                    power_history.clear()                        # clear previous graphs
                    dial_angles[:] = [-90] * NUM_DIALS           # reset dials
                    lever_positions[:] = [LEVER_TRACK_TOP] * NUM_LEVERS  # reset rods to fully out
                continue                                         # skip further click processing for this event
            if power_on and not scram_active:                    # only allow manual rod/dial selection when running
                lever_clicked = False                            # track whether click hit any lever
                for i, x in enumerate(lever_x_positions):        # loop through each lever
                    if (x - 10 < mx < x + 10 and                # within lever handle horizontal range
                        lever_positions[i] < my <               # and vertical range of handle
                        lever_positions[i] + LEVER_HEIGHT):
                        lever_clicked = True                     # mark that a lever was hit
                        lever_active[:] = [False] * NUM_LEVERS   # de-activate all levers
                        dial_active[:]  = [False] * NUM_DIALS    # de-activate all dials
                        selected_dial   = None                   # clear dial selection
                        selected_lever  = i                      # record which lever is now selected
                        lever_active[i] = True                   # mark this lever active
                        lever_activation_time[i] = now           # track activation time for cooldown
                        break                                    # exit lever hit-test loop
                if not lever_clicked:                            # if click did not hit a lever
                    for i, (cx, cy) in enumerate(dial_centers):  # test each dial hit area
                        if math.hypot(mx - cx, my - cy) < DIAL_RADIUS + 10:  # within circle radius + margin
                            lever_active[:] = [False] * NUM_LEVERS           # deactivate levers
                            dial_active[:]  = [False] * NUM_DIALS            # deactivate all dials
                            selected_lever  = None                           # clear lever selection
                            selected_dial   = i                              # set which dial is selected
                            dial_active[i]  = True                           # mark dial active
                            dial_activation_time[i] = now                    # save activation time
                            break                                            # exit dial loop

    keys = pygame.key.get_pressed()                             # snapshot of all keyboard keys this frame
    if power_on and not scram_active:                           # keyboard control only active when running
        if selected_lever is not None and lever_active[selected_lever]:     # if a lever is selected and active
            if keys[pygame.K_UP]:                               # up arrow key pressed
                lever_positions[selected_lever] -= 2            # move lever 2 pixels up (toward fully out)
                lever_activation_time[selected_lever] = now     # refresh last activity time
            elif keys[pygame.K_DOWN]:                           # down arrow key pressed
                lever_positions[selected_lever] += 2            # move lever 2 pixels down (toward fully in)
                lever_activation_time[selected_lever] = now     # refresh activity time
            lever_positions[selected_lever] = max(              # clamp final lever y within allowed track
                LEVER_TRACK_TOP,
                min(LEVER_TRACK_BOTTOM - LEVER_HEIGHT, lever_positions[selected_lever])
            )
        if selected_dial is not None and dial_active[selected_dial]:        # if a dial is selected and active
            if keys[pygame.K_LEFT]:                            # left arrow rotates dial counter-clockwise
                dial_angles[selected_dial] -= 2
                dial_activation_time[selected_dial] = now
            elif keys[pygame.K_RIGHT]:                         # right arrow rotates dial clockwise
                dial_angles[selected_dial] += 2
                dial_activation_time[selected_dial] = now

    for i in range(NUM_LEVERS):                                # iterate through all levers
        if lever_active[i] and (now - lever_activation_time[i]) > COOLDOWN_TIME:
            lever_active[i] = False                            # turn off lever if no activity for COOLDOWN_TIME
            if selected_lever == i:
                selected_lever = None                          # clear selection if that lever timed out

    for i in range(NUM_DIALS):                                 # iterate through all dials
        if dial_active[i] and (now - dial_activation_time[i]) > COOLDOWN_TIME:
            dial_active[i] = False                             # deactivate dial after inactivity
            if selected_dial == i:
                selected_dial = None                           # clear dial selection

    if scram_active:                                           # if SCRAM is in progress
        for i in range(NUM_LEVERS):                            # visually force all rods fully inserted
            lever_positions[i] = LEVER_TRACK_BOTTOM - LEVER_HEIGHT

    if power_on and not scram_active:                          # run reactor physics only when powered and not SCRAM
        rho_fine, rho_coarse = 0.00310, 0.01250                # worth of fine and coarse rods in Δk/k units
        lever_worths = [rho_fine, rho_fine, rho_fine, rho_coarse]  # worth mapping per lever index
        rho_new = 0.0                                          # start with zero total reactivity
        travel = (LEVER_TRACK_BOTTOM - LEVER_TRACK_TOP - LEVER_HEIGHT)  # full travel distance for fraction calc
        for i in range(NUM_LEVERS):                            # accumulate contributions from each rod
            if travel > 0:
                frac_in = (LEVER_TRACK_BOTTOM - LEVER_HEIGHT - lever_positions[i]) / travel  # compute fractional insertion
            else:
                frac_in = 0.0                                  # fallback if travel ever zero
            frac_in = max(0.0, min(1.0, frac_in))              # clamp fractional insertion to [0,1]
            worth_fraction = rod_curve(frac_in)                # nonlinear worth fraction f(x)
            rho_new -= worth_fraction * lever_worths[i]        # inserting rods (increasing frac_in) reduces ρ (more negative)
        RHO_POS_LIMIT, RHO_NEG_LIMIT = 0.95 * beta, -0.2       # hard bounds to keep ρ away from β and too negative
        rho_new = max(RHO_NEG_LIMIT, min(RHO_POS_LIMIT, rho_new))  # clamp ρ within limits

        if abs(rho_new - rho_prev) > 1e-5:                     # if reactivity changed beyond small threshold
            denom_j = (beta - rho_new)                         # denominator for prompt jump factor
            if abs(denom_j) < 1e-6:                            # avoid division by nearly zero near β
                denom_j = 1e-6 if denom_j >= 0 else -1e-6      # enforce tiny nonzero value keeping sign
            raw_jump = (beta - rho_prev) / denom_j             # ideal jump factor n_after/n_before from theory
            raw_jump = max(0.8, min(1.2, raw_jump))            # clamp jump to at most ±20% to keep simulator stable
            jump_factor = 1.0 + 0.5 * (raw_jump - 1.0)         # blend 50% toward no-jump to soften behaviour
            n *= jump_factor                                   # apply instantaneous jump to neutron population
            if n < 1e-12:                                      # enforce tiny positive floor
                n = 1e-12

        substeps = 12                                          # number of RK4 substeps per frame to handle stiffness
        dt_sub = dt / max(1, substeps)                         # small substep time
        for _ in range(substeps):                              # integrate over dt using multiple small steps
            n, C = rk4_step((n, C), dt_sub, rho_new)           # advance (n,C) by dt_sub with RK4
            if not (n == n) or n < 1e-12:                      # check for NaN or negative n
                n = 1e-12
            if not (C == C) or C < 0.0:                        # check for NaN or negative C
                C = 0.0

        current_power = P_SCALE * n                            # map neutron population n to power display P
        current_power = max(0.0, min(MAX_POWER, current_power))# clamp P to [0, MAX_POWER]
        power_history.append(current_power)                    # append newest power sample to history

        if current_power > 4.5 or current_power < 0.5:         # if power leaves safe band [0.5,4.5] W
            scram_active = True                                # trigger SCRAM
            scram_start_time = now                             # record SCRAM time
            selected_lever = selected_dial = None              # clear current selections
            lever_active[:] = [False] * NUM_LEVERS             # deactivate levers
            dial_active[:]  = [False] * NUM_DIALS              # deactivate dials

        rho_prev = rho_new                                     # remember this ρ for next-frame prompt jump
    else:
        n = current_power / P_SCALE                            # if not running, keep n consistent with power
        C = (beta / (Lambda * lam)) * n                        # keep C at steady-state relative to n

    pygame.draw.rect(screen, GREY, (GRAPH_X, GRAPH_Y, GRAPH_WIDTH, GRAPH_HEIGHT))   # draw grey graph background box

    linear_width = GRAPH_WIDTH // 2                            # width for each panel (linear and log)
    right_start  = GRAPH_X + linear_width                      # x where log panel starts

    if len(power_history) > 1:                                 # only draw curves if we have more than one sample
        history_list = list(power_history)                     # copy deque to list for index-based access
        max_points   = len(history_list)                       # number of stored points

        pts = []                                               # list for linear power polyline
        for i, p in enumerate(history_list):                   # iterate through stored power samples
            y = GRAPH_Y + (i / max_points) * GRAPH_HEIGHT      # map index to y (older points higher up)
            x = GRAPH_X + (p / MAX_POWER) * linear_width       # map power 0..MAX_POWER to left panel x range
            pts.append((x, y))                                 # store this point
        pygame.draw.lines(screen, RED, False, smooth_points(pts, steps=8), 2)  # draw smoothed red linear curve

        log_pts = []                                           # list for log power polyline
        for i, p in enumerate(history_list):
            p_clamped = max(p, P_MIN_LOG)                      # enforce minimum value to avoid log10(0)
            log_p = math.log10(p_clamped)                      # compute log10(P)
            frac = (log_p - LOG_MIN) / (LOG_MAX - LOG_MIN)     # normalize log value to [0,1] across log range
            frac = max(0.0, min(1.0, frac))                    # clamp normalized value
            y = GRAPH_Y + (i / max_points) * GRAPH_HEIGHT      # same y mapping as linear
            x = right_start + frac * linear_width              # map normalized log position to right panel x
            log_pts.append((x, y))                             # store point
        pygame.draw.lines(screen, YELLOW, False, log_pts, 2)   # draw yellow log curve (no smoothing needed)

    pygame.draw.rect(screen, BLACK, (GRAPH_X, GRAPH_Y, linear_width, GRAPH_HEIGHT), 2)      # border around linear panel
    pygame.draw.rect(screen, BLACK, (right_start, GRAPH_Y, linear_width, GRAPH_HEIGHT), 2)  # border around log panel

    font = pygame.font.SysFont(None, 26)                       # main label font
    text = font.render(f"Power: {current_power:.2f} W", True, BLACK)  # current power numeric display
    screen.blit(text, (GRAPH_X, GRAPH_Y - 35))                 # draw at top-left above graph

    font_small = pygame.font.SysFont(None, 20)                 # smaller font for panel titles
    linear_label = font_small.render("Linear Power (0–5 W)", True, BLACK)  # label for left panel
    log_label    = font_small.render("Log Power (10^-3 to 10^1 W)", True, BLACK)  # label for right panel
    screen.blit(linear_label,
                linear_label.get_rect(center=(GRAPH_X + linear_width // 2, GRAPH_Y - 15)))  # center above left panel
    screen.blit(log_label,
                log_label.get_rect(center=(right_start + linear_width // 2, GRAPH_Y - 15))) # center above right panel

    y_label = font.render("Time (index)", True, BLACK)         # y-axis label text
    y_label = pygame.transform.rotate(y_label, 90)             # rotate label to vertical orientation
    screen.blit(y_label,
                y_label.get_rect(center=(GRAPH_X - 55, GRAPH_Y + GRAPH_HEIGHT // 2)))  # place beside left side
    draw_ticks_and_labels()                                    # draw axes ticks after background/labels

    pygame.draw.rect(screen, GREEN if power_on else GREY, power_button_rect)  # draw POWER button filled
    label = font.render("POWER", True, BLACK)                  # POWER text
    screen.blit(label, label.get_rect(center=power_button_rect.center))  # center text on button

    scram_rect = pygame.Rect(MARGIN_LEFT // 2,
                             LEVER_TRACK_TOP + (LEVER_TRACK_BOTTOM - LEVER_TRACK_TOP - 60) // 2,
                             60, 60)                            # define square SCRAM button rectangle
    pygame.draw.rect(screen, SCRAM_RED if scram_active else GREY, scram_rect)  # SCRAM button color depends on state
    scram_label = font.render("SCRAM", True, BLACK)            # SCRAM text
    screen.blit(scram_label,
                scram_label.get_rect(center=(scram_rect.centerx,
                                             scram_rect.centery + 60 + 25)))  # place text label under button

    for i in range(NUM_LEVERS):                                # draw each lever assembly
        x, y = lever_x_positions[i], lever_positions[i]        # current lever x and y
        pygame.draw.rect(screen, GREY,
                         (x - 5, LEVER_TRACK_TOP,
                          10, LEVER_TRACK_BOTTOM - LEVER_TRACK_TOP))   # draw vertical grey track behind lever
        color = (LEVER_INACTIVE_COLOR if scram_active or not power_on  # choose handle color by state
                 else (YELLOW if lever_active[i] else BLUE))
        pygame.draw.rect(screen, color, (x - 10, y, 20, LEVER_HEIGHT)) # draw lever handle rectangle
        label = font.render(lever_names[i], True, BLACK)       # lever name text
        screen.blit(label, label.get_rect(center=(x, LEVER_TRACK_BOTTOM + 25)))  # centered below track

    spectrum_x      = lever_x_positions[-1] + 60               # x-position to the right of last lever for reactivity scale
    spectrum_top    = LEVER_TRACK_TOP                          # top y of scale
    spectrum_bottom = LEVER_TRACK_BOTTOM                       # bottom y of scale
    spectrum_height = spectrum_bottom - spectrum_top           # total vertical extent
    minus_label = font.render("-", True, BLACK)                # "-" label at top for negative reactivity
    screen.blit(minus_label,
                minus_label.get_rect(center=(spectrum_x, spectrum_top - 15)))  # place above scale
    plus_label = font.render("+", True, BLACK)                 # "+" label at bottom for positive reactivity
    screen.blit(plus_label,
                plus_label.get_rect(center=(spectrum_x, spectrum_bottom + 15))) # place below scale
    react_label = font.render("Reactivity", True, BLACK)       # vertical label "Reactivity"
    screen.blit(react_label,
                react_label.get_rect(center=(spectrum_x + 40,
                                             spectrum_top + spectrum_height // 2)))  # place to the right of scale

    for i in range(NUM_DIALS):                                 # draw each dial (Lin Zoom, Log Zoom)
        cx, cy = dial_centers[i]                               # center coordinates of dial
        pygame.draw.circle(screen, GREY, (cx, cy), DIAL_RADIUS + 7)  # outer bezel circle
        color = (DIAL_INACTIVE_COLOR if not power_on           # dial face color depends on power and activity
                 else (YELLOW if dial_active[i] else DIAL_INACTIVE_COLOR))
        pygame.draw.circle(screen, color, (cx, cy), DIAL_RADIUS)      # inner circle of dial face
        draw_dial_ticks(cx, cy, DIAL_RADIUS)                   # draw tick marks and labels around dial
        angle_rad = math.radians(dial_angles[i])               # convert dial angle from degrees to radians
        line_x = cx + DIAL_RADIUS * math.cos(angle_rad)        # compute pointer tip x using cos(θ)
        line_y = cy + DIAL_RADIUS * math.sin(angle_rad)        # compute pointer tip y using sin(θ)
        pygame.draw.line(screen, WHITE, (cx, cy), (line_x, line_y), 4) # draw white pointer line
        label = font.render(dial_names[i], True, BLACK)        # dial caption text
        screen.blit(label,
                    label.get_rect(center=(cx, cy - DIAL_RADIUS - 40)))  # position caption above dial

    pygame.display.flip()                                      # swap display buffers, present frame
    clock.tick(60)                                             # limit frame rate to ~60 FPS

pygame.quit()                                                  # uninitialize all pygame modules
sys.exit()                                                     # exit Python interpreter cleanly
