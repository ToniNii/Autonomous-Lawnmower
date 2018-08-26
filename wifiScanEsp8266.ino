 /*
    HOX: NodeMCU 1.0 
  
    This sketch demonstrates how to scan WiFi networks.
    The API is almost the same as with the WiFi Shield library,
    the most obvious difference being the different file you need to include:
*/
#include "ESP8266WiFi.h"
#include <math.h>


void setup() {
  Serial.begin(115200);

  // Set WiFi to station mode and disconnect from an AP if it was previously connected
  WiFi.mode(WIFI_STA);
  WiFi.disconnect();
  delay(100);
  

  Serial.println("Setup done");
}

void loop() {
  Serial.println("scan start");

  // WiFi.scanNetworks will return the number of networks found
  int n = WiFi.scanNetworks();
  double distance=0;
  int rssi=0;
  double to_exp=0;
  Serial.println("scan done");
  if (n == 0) {
    Serial.println("no networks found");
  } else {
    Serial.print(n);
    Serial.println(" networks found");
    for (int i = 0; i < n; ++i) {
      if (WiFi.SSID(i)=="SSID"){ //|| WiFi.SSID(i)=="iPhone (Mirella)" || WiFi.SSID(i)=="Appo")
        Serial.print(WiFi.SSID(i)+ " ,strength: "+ WiFi.RSSI(i) +"dBm");
        rssi=-WiFi.RSSI(i);
        to_exp=(27.55-(20*log10(2412))+rssi)/20;
        distance=100*pow(10,to_exp);
        Serial.print(" ,distance: ");
        Serial.print((float)(distance),1);
        Serial.println("cm");
      }
    }
  }
  Serial.println("");

  // Wait a bit before scanning again
  delay(3000);
}
