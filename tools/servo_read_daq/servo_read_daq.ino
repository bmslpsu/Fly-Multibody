
#include <Servo.h>

#define daq_pin 3
#define servo_pin 9

int daq_signal;
int deg_val;


Servo myservo; // create servo object

void setup() {
  // 
  myservo.attach(servo_pin); // attach servo
  myservo.write(0);
  delay(1);
}

void loop() {
  daq_signal = analogRead(daq_pin);
  deg_val = map(daq_signal, 0, 1023, 0, 180);
  myservo.write(deg_val);
  delay(1);
}
