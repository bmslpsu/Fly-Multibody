
#include <Servo.h>

#define daq_pin A0
#define servo_pin 9

int daq_signal;
int deg_val;

Servo myservo; // create servo object

void setup() {
  // 
  Serial.begin(9600);
  Serial.println("Initialize");
  myservo.attach(servo_pin); // attach servo
  myservo.write(0);
  delay(5);
}

void loop() {
  daq_signal = analogRead(daq_pin);
  deg_val = map(daq_signal, 0, 1023, 0, 92.2904);
  myservo.write(deg_val);
  Serial.println(deg_val);
  delay(1);
}
