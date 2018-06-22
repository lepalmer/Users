""" 
Code to move the linear actuator and rotation table via serial port.

"""

import serial
import time


def command_ser(handler, command):
    """
    This function writes commands to the serial port and changes the encoding of the input.
    
    Parameters
    ----------
    handler: function, string
        Calls the Serial function from the serial module (part of the pyserial library) and specifies a USB port.
        The USB port must be a string.
        An example of a handler would be:
        ser = serial.Serial('/dev/cu.usbserial-A1069JDR', timeout = 60)

    command: string
        The command or phrase being sent to the machine via the serial connection.
        If you want to move the linear actuator x amount of steps, you write:
        command_ser(handler, '@D x')
        The @ sets the distance for all axes.
        Alternatively, if you just wanted to set the distance for axis 2, write:
        command_ser(handler, 'D 0,x,0,0')
        Another useful command is the HOM command, which positions the linear actuator at a preset 'home.'
        You would write:
        command_ser(handler, 'HOM 1111')

    Returns
    -------
    Nothing

    """
    TERMINATOR = b'\r\n'
    ENCODING = 'utf-8'
    UNICODE_HANDLING = 'replace'

    handler.write(command.encode(ENCODING, UNICODE_HANDLING) + TERMINATOR)
    
    return;


def init_cond(handler, scaling, acc=20000, vel=20000):
    """
    Sets the initial conditions for the linear actuator and rotation table.

    Parameters
    ----------
    handler: function, string
        Calls the Serial function from the serial module (part of the pyserial library) and specifies a USB port.
        An example of a handler would be:
        ser = serial.Serial('/dev/cu.usbserial-A1069JDR', timeout = 60)

    scaling: float
        The scaling factor. All acceleration, velocity, and distance values will be multiplied by the scaling factor.

    acc: float
        The acceleration (recommended values between 5000 steps/sec^2 and 30000 steps/sec^2).
        NOTE: There are 5,000 steps per mm.

    vel: float
        The velocity (recommended values between 5000 steps/sec and 30000 steps/sec).
        NOTE: There are 5,000 steps per mm.

    Returns
    -------
    Nothing

    """

    command_ser(handler, '@LH0') # Disables hard limits
    command_ser(handler, '@SCALE %f' % scaling) # Turns on scaling factor
    command_ser(handler, '@AXSDEF 0') # Defines motor as stepper, not servo
    command_ser(handler, 'DRIVE111') # Turns on drives
    command_ser(handler, '@A %f' % acc) # Sets acceleration (multiplied by scaling factor) 
    command_ser(handler, '@V %f' % vel) # Sets velocity (multiplied by scaling factor)
    
    return;

def move_linAct(handler, amount, axis):
    """
    This will move the linear actuator a certain number of mm.
    You can move either axis or both at once.

    Parameters
    ----------

    handler: function, string
        Calls the Serial function from the serial module (part of the pyserial library) and specifies a USB port.
        An example of a handler would be:
        ser = serial.Serial('/dev/cu.usbserial-A1069JDR', timeout = 60)


    amount: float
        The amount term must be in mmm. Axis 1 (the x axis) has a max travel distance of ~162 mm and Axis
        2 (the y axis) has a max travel distance of ~145 mm; therefore, the amount term should not exceed
        these values.
        NOTE: Make sure you remember to specify a positive or negative amount to indicate direction.

    axis: float, string
        This term specifies which axis you want moved. 1 will move Axis 1, 2 will move Axis 2, and 'both' 
        will move both axes.

    Returns
    -------
    Nothing

    """
    
    command_ser(handler, '@D %f' %(amount*5000)) # there are 5000 steps per mm (determined experimentally)
    
    if axis == 1:   
        command_ser(handler, 'GO 1,0,0,0')
    else:
        if axis == 2:
            command_ser(handler, 'GO 0,1,0,0')
        else:
            if axis == 'both':
                command_ser(handler, 'GO 1,1,0,0')
        
    return;

def move_platform(handler, degs):
    """
    Moves the rotation table a desired amount of degrees.

    Parameters
    ----------

    handler: function, string
        Calls the Serial function from the serial module (part of the pyserial library) and specifies a USB port.
        An example of a handler would be:
        ser = serial.Serial('/dev/cu.usbserial-A1069JDR', timeout = 60)

    degs: float
        How many degrees you want the table to move.

    Return
    ------
    Nothing

    """

    dist = int(degs*3125)

    command_ser(handler, 'D 0,0,%f,0' % dist)
    command_ser(handler, 'GO 0,0,1,0')

    return;

def interval_rot(handler, degs, secs, max_angle=90):
    """
    This will move the rotation table a maximum number of degrees, pasuing at a specified intervals for
    a specified amount of time.

    Parameters
    ----------
    handler: function, string
        Calls the Serial function from the serial module (part of the pyserial library) and specifies a USB port.
        An example of a handler would be:
        ser = serial.Serial('/dev/cu.usbserial-A1069JDR', timeout = 60)

    degs: float
        Refers to the number of degrees you want rotated per interval.

    secs: float
        Refers to the the amount of seconds the platform is pausing at each interval.

    max_angle: float
        The maximum angle you want the rotation table to move.

    Returns
    -------
    degs, intervals

    """
        
    dist = int(degs*3125)
    intervals = int(max_angle/degs) + 1
    
    for x in range(0, intervals):
        command_ser(handler, 'D 0,0,%f,0' % dist)
        command_ser(handler, 'GO 0010')
        time.sleep(secs)
         
    return degs, intervals

def stop(handler):
    """
    Kills motion for both the linear actuator and the rotation table.

    Parameters
    ----------
    handler: function, string
        Calls the Serial function from the serial module (part of the pyserial library) and specifies a USB port.
        An example of a handler would be:
        ser = serial.Serial('/dev/cu.usbserial-A1069JDR', timeout = 60)

    Returns
    -------
    Nothing

    """
    
    command_ser(handler, 'K!')
    
    return;
